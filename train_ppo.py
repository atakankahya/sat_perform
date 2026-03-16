"""
PPO training for UPMSat-2 magnetic ADCS
Matches Pérez-Muñoz et al. (2025) Table 4 hyperparameters

Build sat_sim first:
    mkdir build_python && cd build_python
    cmake .. -DPYTHON_EXECUTABLE=$(python3 -c "import sys; print(sys.executable)")
    cmake --build . --target sat_sim --config Release -j4
    cp sat_sim*.so ..   (or .pyd on Windows)

Then: python train_ppo.py
"""
import os
os.add_dll_directory(r"C:\Program Files\JetBrains\CLion 2025.3.2\bin\mingw\bin")
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.distributions import Categorical
import time


import sat_sim

# ── hyperparameters (Table 4) ──
SEED = 42
NUM_ENVS = 8
NUM_EPISODES = 1000
STEPS_PER_EPISODE = 4000       # 4000 steps * 2s = 8000s sim time
CLIP_EPS = 0.2
ENTROPY_COEF = 2e-4
EXPERIENCE_HORIZON = 1024
MINIBATCH_SIZE = 512
DISCOUNT = 0.99
GAE_LAMBDA = 0.95
UPDATE_EPOCHS = 4              # not specified in paper, standard PPO value
ACTOR_LR = 1e-5
CRITIC_LR = 1e-4

OBS_DIM = sat_sim.SatEnv.obs_dim      # 9
N_ACTIONS = sat_sim.SatEnv.n_actions  # 125

IGRF_DIR = "sensors"


# ── separate actor and critic (matching paper architecture) ──

class Actor(nn.Module):
    def __init__(self, obs_dim, n_actions):
        super().__init__()
        # Table 4: Actor hidden layers 400, 256, 200, 200 (ReLU)
        self.net = nn.Sequential(
            nn.Linear(obs_dim, 400), nn.ReLU(),
            nn.Linear(400, 256),     nn.ReLU(),
            nn.Linear(256, 200),     nn.ReLU(),
            nn.Linear(200, 200),     nn.ReLU(),
            nn.Linear(200, n_actions),  # linear output (logits)
        )

    def forward(self, x):
        return self.net(x)

    def get_dist(self, obs):
        logits = self.forward(obs)
        return Categorical(logits=logits)


class Critic(nn.Module):
    def __init__(self, obs_dim):
        super().__init__()
        # Table 4: Critic hidden layers 400, 200, 200 (ReLU)
        self.net = nn.Sequential(
            nn.Linear(obs_dim, 400), nn.ReLU(),
            nn.Linear(400, 200),     nn.ReLU(),
            nn.Linear(200, 200),     nn.ReLU(),
            nn.Linear(200, 1),       # linear output (value)
        )

    def forward(self, x):
        return self.net(x).squeeze(-1)


# ── vectorized environment ──

class VecEnv:
    def __init__(self, num_envs, igrf_dir):
        self.envs = [sat_sim.SatEnv(igrf_dir) for _ in range(num_envs)]
        self.num_envs = num_envs

    def reset(self):
        obs = np.zeros((self.num_envs, OBS_DIM), dtype=np.float32)
        for i, env in enumerate(self.envs):
            obs[i] = env.reset()
        return obs

    def step(self, actions):
        obs = np.zeros((self.num_envs, OBS_DIM), dtype=np.float32)
        rewards = np.zeros(self.num_envs, dtype=np.float32)
        dones = np.zeros(self.num_envs, dtype=np.float32)
        for i, (env, act) in enumerate(zip(self.envs, actions)):
            o, r, d = env.step(int(act))
            rewards[i] = r
            dones[i] = float(d)
            if d:
                obs[i] = env.reset()
            else:
                obs[i] = o
        return obs, rewards, dones


# ── training loop ──

def train():
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    device = torch.device("cpu")

    envs = VecEnv(NUM_ENVS, IGRF_DIR)
    actor = Actor(OBS_DIM, N_ACTIONS).to(device)
    critic = Critic(OBS_DIM).to(device)
    actor_optim = optim.Adam(actor.parameters(), lr=ACTOR_LR, eps=1e-5)
    critic_optim = optim.Adam(critic.parameters(), lr=CRITIC_LR, eps=1e-5)

    # total training updates
    total_steps = NUM_EPISODES * STEPS_PER_EPISODE
    steps_per_collect = EXPERIENCE_HORIZON * NUM_ENVS
    num_updates = total_steps // steps_per_collect

    # rollout storage
    obs_buf = torch.zeros((EXPERIENCE_HORIZON, NUM_ENVS, OBS_DIM), device=device)
    act_buf = torch.zeros((EXPERIENCE_HORIZON, NUM_ENVS), dtype=torch.long, device=device)
    logp_buf = torch.zeros((EXPERIENCE_HORIZON, NUM_ENVS), device=device)
    rew_buf = torch.zeros((EXPERIENCE_HORIZON, NUM_ENVS), device=device)
    done_buf = torch.zeros((EXPERIENCE_HORIZON, NUM_ENVS), device=device)
    val_buf = torch.zeros((EXPERIENCE_HORIZON, NUM_ENVS), device=device)

    next_obs = torch.tensor(envs.reset(), dtype=torch.float32, device=device)
    next_done = torch.zeros(NUM_ENVS, device=device)

    global_step = 0
    start_time = time.time()

    print(f"Training: {num_updates} updates, {total_steps} total steps")
    print(f"Actor: 400-256-200-200 ReLU, lr={ACTOR_LR}")
    print(f"Critic: 400-200-200 ReLU, lr={CRITIC_LR}")
    print(f"Actions: 125 discrete, Obs: 9-dim normalized")
    print("=" * 70)

    for update in range(1, num_updates + 1):
        # ── collect rollout ──
        for step in range(EXPERIENCE_HORIZON):
            obs_buf[step] = next_obs
            done_buf[step] = next_done

            with torch.no_grad():
                dist = actor.get_dist(next_obs)
                action = dist.sample()
                logprob = dist.log_prob(action)
                value = critic(next_obs)

            act_buf[step] = action
            logp_buf[step] = logprob
            val_buf[step] = value

            o, r, d = envs.step(action.cpu().numpy())
            next_obs = torch.tensor(o, dtype=torch.float32, device=device)
            rew_buf[step] = torch.tensor(r, dtype=torch.float32, device=device)
            next_done = torch.tensor(d, dtype=torch.float32, device=device)
            global_step += NUM_ENVS

        # ── GAE ──
        with torch.no_grad():
            next_value = critic(next_obs)

        advantages = torch.zeros_like(rew_buf)
        lastgaelam = 0
        for t in reversed(range(EXPERIENCE_HORIZON)):
            if t == EXPERIENCE_HORIZON - 1:
                nextnonterminal = 1.0 - next_done
                nextvalues = next_value
            else:
                nextnonterminal = 1.0 - done_buf[t + 1]
                nextvalues = val_buf[t + 1]
            delta = rew_buf[t] + DISCOUNT * nextvalues * nextnonterminal - val_buf[t]
            advantages[t] = lastgaelam = delta + DISCOUNT * GAE_LAMBDA * nextnonterminal * lastgaelam

        returns = advantages + val_buf

        # ── flatten ──
        b_obs = obs_buf.reshape(-1, OBS_DIM)
        b_act = act_buf.reshape(-1)
        b_logp = logp_buf.reshape(-1)
        b_adv = advantages.reshape(-1)
        b_ret = returns.reshape(-1)

        batch_size = EXPERIENCE_HORIZON * NUM_ENVS

        # ── PPO update ──
        for epoch in range(UPDATE_EPOCHS):
            indices = torch.randperm(batch_size, device=device)
            for start in range(0, batch_size, MINIBATCH_SIZE):
                end = start + MINIBATCH_SIZE
                if end > batch_size:
                    break
                mb_idx = indices[start:end]

                # actor loss
                dist = actor.get_dist(b_obs[mb_idx])
                newlogprob = dist.log_prob(b_act[mb_idx])
                entropy = dist.entropy()

                mb_adv = b_adv[mb_idx]
                mb_adv = (mb_adv - mb_adv.mean()) / (mb_adv.std() + 1e-8)

                ratio = torch.exp(newlogprob - b_logp[mb_idx])
                pg_loss1 = -mb_adv * ratio
                pg_loss2 = -mb_adv * torch.clamp(ratio, 1 - CLIP_EPS, 1 + CLIP_EPS)
                pg_loss = torch.max(pg_loss1, pg_loss2).mean()
                ent_loss = entropy.mean()

                actor_loss = pg_loss - ENTROPY_COEF * ent_loss

                actor_optim.zero_grad()
                actor_loss.backward()
                nn.utils.clip_grad_norm_(actor.parameters(), 0.5)
                actor_optim.step()

                # critic loss
                newvalue = critic(b_obs[mb_idx])
                v_loss = 0.5 * ((newvalue - b_ret[mb_idx]) ** 2).mean()

                critic_optim.zero_grad()
                v_loss.backward()
                nn.utils.clip_grad_norm_(critic.parameters(), 0.5)
                critic_optim.step()

        # ── logging ──
        if update % 5 == 0:
            elapsed = time.time() - start_time
            sps = global_step / elapsed
            avg_reward = rew_buf.mean().item()
            print(f"Update {update:4d}/{num_updates} | "
                  f"Steps {global_step:8d} | "
                  f"SPS {sps:.0f} | "
                  f"Avg Reward {avg_reward:.4f} | "
                  f"Time {elapsed:.0f}s")

        # ── checkpoint ──
        if update % 50 == 0:
            os.makedirs("checkpoints", exist_ok=True)
            torch.save({
                'actor': actor.state_dict(),
                'critic': critic.state_dict(),
                'update': update,
                'global_step': global_step,
            }, f"checkpoints/ppo_update_{update}.pt")

            # TorchScript for C++ inference
            scripted = torch.jit.script(actor)
            scripted.save(f"checkpoints/actor_update_{update}_scripted.pt")
            print(f"  -> Checkpoint saved at update {update}")

    # final save
    os.makedirs("checkpoints", exist_ok=True)
    torch.save({
        'actor': actor.state_dict(),
        'critic': critic.state_dict(),
    }, "checkpoints/ppo_final.pt")
    scripted = torch.jit.script(actor)
    scripted.save("checkpoints/actor_final_scripted.pt")
    print("Training complete.")


if __name__ == "__main__":
    train()