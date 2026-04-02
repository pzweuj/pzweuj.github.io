---
title: 构建系统镜像及硬盘快照
tags: default
---

我**对成本极其敏感**，需要**高可用生信环境**的云服务器进行系统镜像和硬盘快照构建，提升云服务器的启动速度。

---

## 核心逻辑：三位一体架构
1.  **系统镜像（Custom Image）：** 存放驱动和 Docker 环境（约 20GB）。
2.  **数据快照（Data Snapshot）：** 存放数据库（hg38等）和 Docker 镜像数据（约 200GB）。
3.  **计算实例（Spot Instance）：** 临时启动，挂载快照，跑完即销毁。

---

## 第一阶段：炼制“系统镜像”（GPU 环境）
**目标：** 制作一个预装好驱动和 Docker 的“模具”。
**成本估算：** 使用按量计费 T4 实例，约 1 小时，耗费约 **3-5 元**。

1.  **开机：** 在腾讯云创建一台 **按量计费** 的 GPU 实例（如 GN7-T4）。镜像选官方 **Ubuntu 22.04 Server**（纯净版）。

2.  **安装基础组件：**
    ```bash
    sudo apt-get update
    sudo apt-get install -y docker.io nfs-common
    ```
3.  **安装 NVIDIA 驱动与 Toolkit：**
    * 使用腾讯云脚本或手动安装驱动。
    * 安装连接工具（核心）：
        ```bash
        # 添加源
        curl -s -L https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg
        # 安装 toolkit
        sudo apt-get install -y nvidia-container-toolkit
        ```
4.  **【关键】配置 Docker 路径重定向：**
    编辑 `/etc/docker/daemon.json`（如果没有则新建）：
    ```json
    {
      "data-root": "/mnt/data/docker_root",
      "default-runtime": "nvidia",
      "runtimes": {
        "nvidia": {
          "path": "nvidia-container-runtime",
          "runtimeArgs": []
        }
      }
    }
    ```
    *注意：此时重启 Docker 会报错，因为 `/mnt/data` 还没挂载，没关系。*
5.  **保存镜像：** 在腾讯云控制台点击“制作镜像”，命名为 `WES_Base_GPU_Image`。
6.  **销毁：** 立即销毁该 GPU 实例。

---

## 第二阶段：填充“数据盘快照”（CPU 环境）
**目标：** 在便宜的机器上把 TB 级的数据和镜像塞进硬盘。
**成本估算：** 按量计费 S5 实例，约 2 小时，耗费约 **1 元**。

1.  **开机：** 创建一台最便宜的 **CPU 实例**，同时购买一块 **200GB 的按量计费云硬盘**。
2.  **挂载硬盘：**
    ```bash
    mkfs.ext4 /dev/vdb  # 格式化新盘
    mkdir -p /mnt/data
    mount /dev/vdb /mnt/data
    mkdir -p /mnt/data/docker_root
    ```
3.  **同步 Docker 配置：**
    重复第一阶段第 4 步的 `daemon.json` 配置，然后重启 Docker：
    ```bash
    sudo systemctl restart docker
    ```
4.  **填充资产：**
    * **拉取生信镜像：** `docker pull google/deepvariant:1.6.0-gpu`（镜像现在物理存储在 `/mnt/data/docker_root`）。
    * **下载数据库：** 在 `/mnt/data/database/` 下放入 `hg38.fa` 等。
5.  **保存快照：** 在“云硬盘”页面，点击该 200GB 硬盘的“**创建快照**”，命名为 `WES_Asset_Snapshot`。
6.  **销毁：** 销毁该 CPU 实例和 200GB 云硬盘（快照已存）。



---

## 第三阶段：正式生产（批量并行）
**目标：** 1 小时内拉起 20 台机器处理 20 个样本。
**成本估算：** 竞价实例，每个样本约 **2-4 元**。

### 1. 自动化启动脚本（UserData）
在批量创建实例的“高级设置”中，输入以下脚本：

```bash
#!/bin/bash
# 1. 自动挂载克隆出来的工具盘
mkdir -p /mnt/data
mount /dev/vdb1 /mnt/data

# 2. 准备临时计算空间（系统盘剩余空间或另挂硬盘）
mkdir -p /mnt/workspace

# 3. 启动 Docker（它会自动寻找 /mnt/data 里的镜像）
systemctl restart docker

# 4. 获取任务并执行（假设你把样本名存在了实例名称里）
SAMPLE_NAME=$(curl -s http://metadata.tencentyun.com/latest/meta-data/instance-name)

# 运行 Docker 分析
docker run --rm --gpus all \
  -v /mnt/data/database:/db \
  -v /mnt/workspace:/work \
  google/deepvariant:1.6.0-gpu \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WES \
  --ref=/db/hg38.fa \
  --reads=/work/$SAMPLE_NAME.bam \
  --output_vcf=/work/$SAMPLE_NAME.vcf

# 5. 上传结果并自杀
coscmd upload /mnt/workspace/$SAMPLE_NAME.vcf /results/
shutdown -h now
```

### 2. 控制台操作要点
* **付费模式：** 必须选“**竞价实例**”。
* **镜像：** 选你做的 `WES_Base_GPU_Image`。
* **数据盘：** 点击“添加云硬盘”，选择“**从快照创建**”，选中 `WES_Asset_Snapshot`。
* **自动释放：** 勾选“**关机后自动销毁**”。

---

## 方案总结：如何确保“资产在快照里”？

1.  **路径验证：** 只要你的 `daemon.json` 里的 `data-root` 指向了 `/mnt/data`，那么 `docker pull` 的每一比特数据都绝对在数据盘里。
2.  **只读安全：** 在生产阶段，你可以用 `mount -o ro /dev/vdb1 /mnt/data` 挂载，这样 20 台机器不管怎么跑，都绝对无法破坏你的“资产快照”。
3.  **结果导向：** 记住，**系统盘是消耗品，快照是固定资产，COS 是仓库**。

