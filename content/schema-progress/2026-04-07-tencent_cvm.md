---
title: 腾讯云镜像构建与数据传输指南
tags: default
---

## 背景

在进行基因组数据分析任务时，我需要频繁调起云服务器执行比对和变异检测。为了优化成本，采用竞价实例配合自定义镜像的方案，实现快速启动和自动化部署。

## 镜像构建思路

每个任务调起一台竞价计费的云服务器执行。优先申请A10或T4 GPU实例以加速分析任务，但由于AI热潮导致GPU资源紧张，需要准备CPU实例作为备选方案。

构建自定义镜像的核心流程：

1. 安装NVIDIA驱动（GPU实例）
2. 安装Docker并迁移数据目录到数据盘
3. 配置nvidia-container-runtime
4. 拉取Parabricks镜像并测试

### 1. 申请基础实例

申请一台GPU实例进行镜像构建，同时添加一块300GB**高性能云硬盘**存储数据库和Docker镜像。

### 2. NVIDIA驱动安装

腾讯云的GPU实例在申请时，可以勾选自动安装驱动。进入服务器后会提示驱动正在安装，在完成安装后，使用下面命令验证即可。

```bash
# 验证安装
nvidia-smi
```

### 3. Docker安装与配置

```bash
# 安装Docker
curl -fsSL https://get.docker.com | sh

# 将当前用户加入docker组
sudo usermod -aG docker $USER

# 重新登录后验证
docker --version
```

### 4. 云硬盘挂载

```bash
# 查看磁盘信息
lsblk

# 假设数据盘为 /dev/vdb，进行分区和格式化
sudo fdisk /dev/vdb
# 输入: n -> p -> 1 -> 回车 -> 回车 -> w

# 格式化为ext4
sudo mkfs.ext4 /dev/vdb1

# 设置卷标（重要：便于自动化挂载）
sudo e2label /dev/vdb1 datadisk

# 创建挂载点
sudo mkdir -p /data

# 挂载
sudo mount /dev/disk/by-label/datadisk /data

# 设置权限
sudo chown -R $USER:$USER /data

# 配置开机自动挂载
echo "LABEL=datadisk /data ext4 defaults 0 2" | sudo tee -a /etc/fstab
```

### 5. Docker数据目录迁移

```bash
# 停止Docker服务
sudo systemctl stop docker
sudo systemctl stop docker.socket

# 迁移数据目录
sudo rsync -aP /var/lib/docker/ /data/docker/

# 修改Docker配置
sudo tee /etc/docker/daemon.json <<EOF
{
  "data-root": "/data/docker",
  "storage-driver": "overlay2"
}
EOF

# 重启Docker
sudo systemctl start docker

# 验证新路径
docker info | grep "Docker Root Dir"
```

### 6. 安装nvidia-container-runtime

```bash
# 添加NVIDIA Container Toolkit源
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg

curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
  sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
  sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list

# 安装
sudo apt-get update
sudo apt-get install -y nvidia-container-toolkit

# 配置Docker运行时
sudo nvidia-ctk runtime configure --runtime=docker

# 重启Docker
sudo systemctl restart docker

# 验证GPU可用
docker run --rm --gpus all nvidia/cuda:12.0-base-ubuntu22.04 nvidia-smi
```

### 7. Parabricks镜像拉取与测试

```bash
# 拉取镜像
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1

# 测试
docker run --rm --gpus all nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1 pbrun version
```

### 8. 创建镜像和快照

在腾讯云控制台执行：

1. **创建自定义镜像**：仅备份系统盘
2. **创建云硬盘快照**：备份300GB数据盘

完成后销毁构建实例。

## 数据库传输方案

数据库已公开托管于Huggingface。由于CoLab下载速度不理想，采用以下传输链路：

```
hf-mirror → 家宽下载 → 腾讯云COS → 云服务器
```

**成本优化**：

- 入流量免费（家宽下载、上传至COS）
- 内网传输免费（从COS下载到云服务器。为什么不直接直接下载到云硬盘：因为项目还在构建阶段，没有必要这么早交额外的储存费）
- 增强型SSD云硬盘快照：约40元/月/300GB
- 数据持久化存储于快照，任务启动时从快照创建新盘

### 数据下载到COS

```bash
# 安装huggingface-cli
pip install huggingface_hub

# 使用镜像站下载
HF_ENDPOINT=https://hf-mirror.com huggingface-cli download \
  your-org/your-database \
  --local-dir ./database \
  --local-dir-use-symlinks False

# 安装coscmd并上传
pip install coscmd
coscmd config -a <SecretId> -s <SecretKey> -b <BucketName> -r ap-guangzhou
coscmd upload -r ./database /database/
```

### 任务启动时恢复数据

```bash
# 从快照创建的云硬盘自动挂载到 /data
# 如需手动挂载
sudo mount /dev/disk/by-label/datadisk /data

# 数据已在 /data 目录就绪
```

## 实例模板配置

在腾讯云控制台创建实例模板：

1. 选择自定义镜像
2. 选择实例规格（GPU优先，CPU备选）
3. 添加数据盘快照
4. 配置网络和安全组
5. 设置竞价实例上限价格

后续可通过API或控制台一键部署，实现任务自动化。

---

> 💡 **提示**：增强型SSD相比高性能云硬盘IOPS提升约2倍，适合生信分析的高IO场景。

