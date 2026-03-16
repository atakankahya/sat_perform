"""
Evaluate trained PPO actor: save CSV + plot angular velocity, PWM, reward
Usage: python eval_agent.py [checkpoint_path]
"""
import sys
import os
if sys.platform == "win32":
    mingw_path = r"C:\Program Files\JetBrains\CLion 2025.3.2\bin\mingw\bin"
    if os.path.exists(mingw_path):
        os.add_dll_directory(mingw_path)

import sys
import numpy as np
import torch
import torch.nn as nn
import csv
import matplotlib.pyplot as plt

import sat_sim


class Actor(nn.Module):
    def __init__(self, obs_dim, n_actions):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(obs_dim, 400), nn.ReLU(),
            nn.Linear(400, 256),     nn.ReLU(),
            nn.Linear(256, 200),     nn.ReLU(),
            nn.Linear(200, 200),     nn.ReLU(),
            nn.Linear(200, n_actions),
        )

    def forward(self, x):
        return self.net(x)


OBS_DIM = sat_sim.SatEnv.obs_dim
N_ACTIONS = sat_sim.SatEnv.n_actions

PWM_LEVELS = [-500.0, -200.0, 0.0, 200.0, 500.0]

# denormalization bounds (Table 3)
WERR_MAX = 0.25


def decode_action(action):
    a2 = action % 5
    a1 = (action // 5) % 5
    a0 = action // 25
    return PWM_LEVELS[a0], PWM_LEVELS[a1], PWM_LEVELS[a2]


def evaluate(checkpoint_path, output_csv="eval_results.csv"):
    actor = Actor(OBS_DIM, N_ACTIONS)
    ckpt = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    if "actor" in ckpt:
        actor.load_state_dict(ckpt["actor"])
    else:
        actor.load_state_dict(ckpt)
    actor.eval()

    env = sat_sim.SatEnv("sensors")
    obs = env.reset()

    times = []
    werr_x, werr_y, werr_z = [], [], []
    pwm_xs, pwm_ys, pwm_zs = [], [], []
    rewards = []
    actions_list = []

    step = 0
    total_reward = 0.0
    done = False

    while not done:
        obs_tensor = torch.tensor(obs, dtype=torch.float32).unsqueeze(0)
        with torch.no_grad():
            logits = actor(obs_tensor)
            action = logits.argmax(dim=1).item()

        pwm_x, pwm_y, pwm_z = decode_action(action)
        obs_next, reward, done = env.step(action)
        total_reward += reward

        t = step * 2.0
        times.append(t)

        # obs[6,7,8] are normalized omega errors in [0,1]
        # denormalize: werr = obs * WERR_MAX
        # but these are |w_desired - w|, we want actual omega
        # w = w_desired - werr  (approximately)
        # w_desired = [0, 0, 0.1]
        werr_x.append(obs[6] * WERR_MAX)
        werr_y.append(obs[7] * WERR_MAX)
        werr_z.append(obs[8] * WERR_MAX)

        pwm_xs.append(pwm_x)
        pwm_ys.append(pwm_y)
        pwm_zs.append(pwm_z)
        rewards.append(reward)
        actions_list.append(action)

        obs = obs_next
        step += 1

    # ── save CSV ──
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s", "werr_x", "werr_y", "werr_z",
                         "pwm_x_ms", "pwm_y_ms", "pwm_z_ms",
                         "reward", "action"])
        for i in range(len(times)):
            writer.writerow([times[i], werr_x[i], werr_y[i], werr_z[i],
                             pwm_xs[i], pwm_ys[i], pwm_zs[i],
                             rewards[i], actions_list[i]])

    # ── print summary ──
    print(f"Episode: {step} steps, {step*2.0:.0f}s sim time")
    print(f"Total reward: {total_reward:.4f}")
    print(f"Final |omega_error|: x={werr_x[-1]:.4f} y={werr_y[-1]:.4f} z={werr_z[-1]:.4f} rad/s")
    print(f"CSV saved to {output_csv}")

    unique, counts = np.unique(actions_list, return_counts=True)
    print(f"\nTop 10 actions:")
    top_idx = np.argsort(-counts)[:10]
    for i in top_idx:
        px, py, pz = decode_action(unique[i])
        print(f"  Action {unique[i]:3d} ({px:+.0f},{py:+.0f},{pz:+.0f})ms "
              f"-> {counts[i]}x ({100*counts[i]/len(actions_list):.1f}%)")

    # ── plots ──
    times = np.array(times)

    fig, axes = plt.subplots(4, 1, figsize=(12, 14), sharex=True)

    # 1. Angular velocity error per axis
    ax = axes[0]
    ax.plot(times, werr_x, label="|ωx_err|", linewidth=0.8)
    ax.plot(times, werr_y, label="|ωy_err|", linewidth=0.8)
    ax.plot(times, werr_z, label="|ωz_err|", linewidth=0.8)
    ax.axhline(y=0, color="k", linestyle="--", alpha=0.3)
    ax.set_ylabel("Angular velocity error [rad/s]")
    ax.set_title("Angular Velocity Error (|ω_desired - ω|)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 2. PWM commands per axis
    ax = axes[1]
    ax.step(times, pwm_xs, label="PWM_X", linewidth=0.6, alpha=0.8)
    ax.step(times, pwm_ys, label="PWM_Y", linewidth=0.6, alpha=0.8)
    ax.step(times, pwm_zs, label="PWM_Z", linewidth=0.6, alpha=0.8)
    ax.set_ylabel("PWM duty cycle [ms]")
    ax.set_title("Magnetorquer PWM Commands")
    ax.set_ylim(-600, 600)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 3. Reward over time
    ax = axes[2]
    ax.plot(times, rewards, linewidth=0.8, color="green")
    ax.set_ylabel("Filtered Reward")
    ax.set_title("Reward per Step")
    ax.grid(True, alpha=0.3)

    # 4. Action distribution histogram
    ax = axes[3]
    ax.hist(actions_list, bins=125, range=(0, 125), color="steelblue", edgecolor="none")
    ax.set_xlabel("Action index (0-124)")
    ax.set_ylabel("Count")
    ax.set_title("Action Distribution")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("eval_plots.png", dpi=150)
    print(f"Plots saved to eval_plots.png")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        ckpt = sys.argv[1]
    else:
        ckpt = "checkpoints/ppo_final.pt"

    if not os.path.exists(ckpt):
        print(f"Checkpoint not found: {ckpt}")
        if os.path.exists("checkpoints"):
            print("Available:")
            for f in sorted(os.listdir("checkpoints")):
                print(f"  checkpoints/{f}")
        sys.exit(1)

    evaluate(ckpt)