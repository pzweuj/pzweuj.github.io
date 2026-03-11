---
title: 同步Docker镜像到阿里云 ACR加速镜像拉取
tags: coding
---

> 在国内服务器使用 Docker 镜像时，网络问题一直是痛点。本文介绍一种方案：通过 GitHub Actions 将 gcr.io、ghcr.io、quay.io、Docker Hub 等国外镜像源的镜像同步到阿里云 ACR，实现稳定高速的镜像拉取。

## 方案概述

### 核心思路

1. **阿里云 ACR 免费版**：申请阿里云容器镜像服务个人版，可创建 3 个命名空间和 300 个镜像仓库
2. **GitHub Actions**：利用 GitHub 免费的 CI/CD 能力，自动拉取国外镜像并推送到 ACR
3. **images.txt 配置**：只需在文件中列出需要同步的镜像，即可触发自动同步
4. **自定义域名（可选）**：通过 Cloudflare Worker 部署反向代理，简化超长的 ACR 地址

### 方案优势

- **完全免费**：阿里云 ACR 个人版免费，GitHub Actions 每月有 2000 分钟免费额度
- **自动化**：修改 images.txt 即可自动触发同步
- **增量同步**：已存在的镜像会自动跳过，避免重复操作
- **重试机制**：对网络不稳定的镜像源有专门的重试策略

---

## 一、阿里云 ACR 配置

### 1.1 创建 ACR 实例

1. 访问 [阿里云容器镜像服务](https://cr.console.aliyun.com/)
2. 点击「开通服务」，选择「个人版」
3. 填写相关信息，完成开通

### 1.2 创建命名空间

1. 在 ACR 控制台左侧点击「命名空间」
2. 点击「创建命名空间」，输入名称（如 `mirror`）
3. 创建「创建镜像仓库」，注意设定为“公开”，便于拉取；或保持设定为“私人”，拉取时先进行`docker login`

### 1.3 获取 ACR 访问凭证

1. 进入 [阿里云容器镜像服务控制台](https://cr.console.aliyun.com/)
2. 在左侧点击「访问凭证」
3. 点击「设置固定密码」，设置一个用于登录 ACR 的密码
4. 保存好这个密码，以及你的阿里云账号 ID（通常是你的阿里云登录邮箱或手机号）

---

## 二、GitHub 仓库配置

### 2.1 创建私有仓库

1. 登录 GitHub，创建新仓库（如 `docker-mirror`）
2. 选择 **Private** 私有仓库
3. 不需要初始化 README

### 2.2 创建文件结构

在仓库中创建以下文件：

```
docker-mirror/
├── .github/
│   └── workflows/
│       └── sync_acr.yml
├── sync.sh
└── images.txt
```

### 2.3 配置 GitHub Secrets

进入仓库设置 → Secrets and variables → Actions，添加以下 secrets：

| Secret 名称 | 值 |
|------------|---|
| `ACR_REGISTRY` | 例如 `xxxxxxxxxxx.cn-guangzhou.personal.cr.aliyuncs.com`（根据你的地域选择） |
| `ACR_USERNAME` | 阿里云账号 ID（登录邮箱或手机号） |
| `ACR_PASSWORD` | ACR 固定密码 |
| `ACR_NAMESPACE` | 你创建的命名空间名称 |

### 2.4 创建 images.txt

在仓库根目录创建 `images.txt` 文件，每行一个镜像：

```
ghcr.io/dexidp/dex:latest
quay.io/prometheus/prometheus:v2.45.0
nginx:latest
redis:7-alpine
```

**说明**：
- 以 `#` 开头的行会被忽略（注释）
- 支持 ghcr.io、quay.io、Docker Hub 等多种镜像源
- 不带标签的镜像默认使用 `:latest`

---

## 三、核心代码

### 3.1 GitHub Actions Workflow

创建 `.github/workflows/sync_acr.yml`：

```yaml
name: Sync Docker Images to Aliyun ACR
on:
  push:
    branches:
      - main
    paths:
      - 'images.txt'
      - '.github/workflows/sync_acr.yml'
  workflow_dispatch:

jobs:
  sync-images:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout repository
        uses: actions/checkout@v4

      - name: Login to Aliyun Container Registry
        uses: docker/login-action@v3
        with:
          registry: ${{ secrets.ACR_REGISTRY }}
          username: ${{ secrets.ACR_USERNAME }}
          password: ${{ secrets.ACR_PASSWORD }}

      - name: Execute sync script
        run: bash "${{ github.workspace }}/sync.sh"
        env:
          ACR_REGISTRY: ${{ secrets.ACR_REGISTRY }}
          ACR_NAMESPACE: ${{ secrets.ACR_NAMESPACE }}
```

### 3.2 同步脚本

创建 `sync.sh`（完整代码见仓库），核心功能包括：
- 读取 images.txt 中的镜像列表
- 自动识别镜像来源（ghcr.io/quay.io/Docker Hub）
- 检测已存在的镜像并跳过
- 失败重试机制（ghcr.io 和 quay.io 最多重试 4 次）
- 推送完成后自动清理本地镜像

```bash
#!/bin/bash
set -eu # -u: 遇到未定义变量报错；-e: 任何命令失败立即退出
 
IMAGES_FILE="images.txt"
SUCCESS_COUNT=0
FAILED_COUNT=0
SKIPPED_COUNT=0
FAILED_IMAGES=()

# 重试函数
retry_command() {
    local max_attempts=2
    local delay=5
    local attempt=1
    local command="$@"
    
    # 对于特殊镜像仓库，调整重试策略
    if [[ "$command" == *"ghcr.io"* ]] || [[ "$command" == *"quay.io"* ]]; then
        max_attempts=4
        delay=8
    fi
    
    while [ $attempt -le $max_attempts ]; do
        echo "🔄 Attempt $attempt/$max_attempts"
        if eval "$command"; then
            echo "✅ Command succeeded"
            return 0
        else
            local exit_code=$?
            echo "❌ Command failed with exit code: $exit_code"
            
            if [ $attempt -lt $max_attempts ]; then
                echo "⏳ Retrying in ${delay} seconds..."
                sleep $delay
                # 对于网络相关的镜像仓库，使用递增延迟
                if [[ "$command" == *"ghcr.io"* ]] || [[ "$command" == *"quay.io"* ]]; then
                    delay=$((delay + 3))
                else
                    delay=$((delay * 2))
                fi
            fi
            attempt=$((attempt + 1))
        fi
    done
    
    echo "💥 Command failed after $max_attempts attempts"
    return 1
}
 
# 检查配置文件和镜像列表文件是否存在
if [ ! -f "$IMAGES_FILE" ]; then
    echo "Error: images.txt not found! Please create it with a list of images to sync."
    exit 1
fi
 
 
# 检查必要的配置变量是否已设置
if [ -z "$ACR_REGISTRY" ] || [ -z "$ACR_NAMESPACE" ]; then
    echo "Error: ACR_REGISTRY or ACR_NAMESPACE not set in github variables. Please check your config."
    exit 1
fi

# 检查镜像类型并显示提示信息
check_image_source() {
    local image="$1"
    if [[ "$image" == ghcr.io/* ]]; then
        echo "🐙 GitHub Container Registry image detected: ${image}"
    elif [[ "$image" == quay.io/* ]]; then
        echo "🔴 Quay.io image detected: ${image}"
    elif [[ "$image" == */* ]]; then
        echo "� Dockerg Hub image detected: ${image}"
    else
        echo "🐳 Docker Hub official image detected: ${image}"
    fi
}

# 处理镜像名称重命名，去掉所有前缀，只保留镜像名
process_image_name() {
    local image="$1"
    local processed_image="$image"

    # 去掉 ghcr.io/ 前缀
    if [[ "$image" == ghcr.io/* ]]; then
        processed_image="${image#ghcr.io/}"
    fi

    # 去掉 quay.io/ 前缀
    if [[ "$image" == quay.io/* ]]; then
        processed_image="${image#quay.io/}"
    fi

    # 去掉组织名/用户名前缀，只保留最终镜像名
    # 例如: pzweuj/mapping -> mapping
    if [[ "$processed_image" == */* ]]; then
        processed_image="${processed_image##*/}"
        echo "🔄 Removed registry/org prefix, final name: ${processed_image}" >&2
    fi

    echo "$processed_image"
}
 
echo "Starting Docker image synchronization to ACR..."
echo "Target Registry: ${ACR_REGISTRY}"
echo "Target Namespace: ${ACR_NAMESPACE}"
echo "-----------------------------------"
 
# 遍历 images.txt，逐行处理镜像
while IFS= read -r image; do
    # 跳过空行或以 # 开头的注释行
    if [[ -z "$image" || "$image" =~ ^# ]]; then
        continue
    fi
 
    echo "--- Processing image: ${image} ---"
    
    # 显示镜像来源信息
    check_image_source "$image"
 
    # 分离原始镜像的仓库名和标签
    # 例如：nginx:latest -> original_repo=nginx, original_tag=latest
    # 例如：jenkins/jenkins:lts -> original_repo=jenkins/jenkins, original_tag=lts
    # 处理没有标签的情况，默认使用 latest
    if [[ "$image" == *":"* ]]; then
        original_repo=$(echo "$image" | cut -d ':' -f1)
        original_tag=$(echo "$image" | cut -d ':' -f2)
    else
        original_repo="$image"
        original_tag="latest"
        image="${image}:latest"  # 更新完整镜像名
    fi
    
    # 处理镜像名称重命名（去掉 ghcr.io 和 quay.io 前缀）
    processed_repo=$(process_image_name "$original_repo")
 
    # 构造目标 ACR 完整镜像路径（使用处理后的镜像名）
    target_full_image_path="${ACR_REGISTRY}/${ACR_NAMESPACE}/${processed_repo}:${original_tag}"
 
    echo "Original image full path: ${image}"
    echo "Target ACR image full path: ${target_full_image_path}"
 
    # 检查阿里云仓库是否已有该tag
    # docker manifest inspect 命令用于检查远程 Registry 中的镜像是否存在。
    # 如果已经存在，则跳过本次同步，避免重复操作和不必要的流量消耗。
    if docker manifest inspect "${target_full_image_path}" > /dev/null 2>&1; then
        echo "${target_full_image_path} 已存在于 ACR，跳过本次同步。"
        SKIPPED_COUNT=$((SKIPPED_COUNT + 1))
        echo "-----------------------------------"
        continue # 跳过当前循环的后续步骤
    fi
 
    echo "Image ${target_full_image_path} not found in ACR. Proceeding with sync..."
 
    # 使用错误处理来捕获同步过程中的失败
    if (
        # 拉取原始镜像
        echo "📥 Pulling original image: ${image}..." &&
        retry_command "docker pull \"${image}\"" &&
        
        # 打上阿里云 ACR 的标签
        echo "Tagging image ${image} to ${target_full_image_path}..." &&
        docker tag "${image}" "${target_full_image_path}" &&
        
        # 推送到阿里云 ACR
        echo "Pushing image ${target_full_image_path} to ACR..." &&
        retry_command "docker push \"${target_full_image_path}\""
    ); then
        echo "Successfully synced: ${image} to ${target_full_image_path}"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "Failed to sync: ${image}"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        FAILED_IMAGES+=("${image}")
    fi
    
    # 清理本地拉取和打标签的镜像，释放 GitHub Actions Runner 的磁盘空间
    echo "Cleaning up local images..."
    # 使用 || true 即使删除失败也不会中断脚本，确保后续镜像能继续处理
    docker rmi "${image}" || true
    docker rmi "${target_full_image_path}" || true
    
    echo "-----------------------------------"
 
done < "$IMAGES_FILE"
 
echo "All specified images processed successfully."
echo "Synchronization process finished."
echo ""
echo "=== SYNC SUMMARY ==="
echo "✅ Successfully synced: ${SUCCESS_COUNT} images"
echo "⏭️  Skipped (already exists): ${SKIPPED_COUNT} images"
echo "❌ Failed: ${FAILED_COUNT} images"

if [ ${FAILED_COUNT} -gt 0 ]; then
    echo ""
    echo "Failed images:"
    for failed_image in "${FAILED_IMAGES[@]}"; do
        echo "  - ${failed_image}"
    done
    echo ""
    echo "⚠️  Some images failed to sync. Please check the logs above for details."
    exit 1
else
    echo ""
    echo "🎉 All images processed successfully!"
fi
```

---

## 四、使用方法

### 4.1 触发同步

**方式一：修改 images.txt**

在 images.txt 中添加新镜像，然后 commit 并 push，GitHub Actions 会自动触发。

**方式二：手动触发**

1. 进入仓库的 Actions 页面
2. 点击「Sync Docker Images to Aliyun ACR」
3. 点击「Run workflow」手动执行

### 4.2 查看同步状态

在 Actions 页面可以查看：
- 同步进度和日志
- 成功/失败/跳过的镜像数量
- 详细的错误信息

### 4.3 使用同步后的镜像

同步完成后，使用以下命令拉取：

```bash
# 完整地址
docker pull xxxxxxxxxxx.cn-guangzhou.personal.cr.aliyuncs.com/mirror/nginx:latest

# 如果配置了自定义域名（见下文）
docker pull my-mirror.io/nginx:latest
```

---

## 五、（可选）Cloudflare Worker 加速

如果觉得 ACR 的域名太长，可以配合 Cloudflare Worker 部署反向代理。

### 5.1 优势

- 简短的自定义域名
- 复用 Cloudflare 全球 CDN 加速
- 支持 Docker pull 协议

### 5.2 部署步骤

1. 在 Cloudflare 控制台创建 Worker
2. 部署以下代码（修改 `UPSTREAM_URL` 为你的 ACR 地址）：

```javascript

const UPSTREAM_URL = '阿里云ACR实例URL';
const CUSTOM_DOMAIN = '自定义域名';

addEventListener('fetch', event => {
  event.respondWith(handleRequest(event.request));
});

async function handleRequest(request) {
  const url = new URL(request.url);
  
  // 1. 替换目标 Host 为阿里云 ACR
  url.hostname = UPSTREAM_URL;

  // 2. 构建新请求头，伪装 Host
  const newRequestHeaders = new Headers(request.headers);
  newRequestHeaders.set('Host', UPSTREAM_URL);

  // 3. 构建请求对象
  // ⚠️ CRITICAL CHANGE: redirect: 'manual'
  // 这告诉 Worker 不要自动跟随重定向，而是把 307 响应原样返回给 Podman
  // 这样 Podman 就会直接去连接阿里云 OSS 下载大文件，不消耗 Worker 资源
  const newRequest = new Request(url.toString(), {
    method: request.method,
    headers: newRequestHeaders,
    body: request.body,
    redirect: 'manual' 
  });

  try {
    const response = await fetch(newRequest);
    const responseHeaders = new Headers(response.headers);

    // 4. 确保 Docker API 版本头存在（Podman对此很敏感）
    if (!responseHeaders.has('Docker-Distribution-Api-Version')) {
        responseHeaders.set('Docker-Distribution-Api-Version', 'registry/2.0');
    }

    // 5. 返回响应
    // 如果是 307 Redirect，Podman 会收到 Location 头并自动处理
    return new Response(response.body, {
      status: response.status,
      statusText: response.statusText,
      headers: responseHeaders
    });

  } catch (e) {
    return new Response(JSON.stringify({ error: e.message }), { status: 500 });
  }
}
```

3. 绑定自定义域名（如 `mirror.yourdomain.com`）

---

## 六、进阶策略：使用 GitHub API 自动化管理镜像列表

### 6.1 为什么需要 API 封装

手动编辑 images.txt 虽然简单，但在以下场景可能不够方便：
- 想通过脚本批量添加/删除镜像
- 想从其他镜像源自动同步可用镜像列表
- 想在内部系统中集成镜像管理功能

### 6.2 核心思路

利用 GitHub REST API 来：
1. **读取** 当前 images.txt 内容
2. **修改** 镜像列表（添加/删除）
3. **提交** 更改到仓库（自动触发 GitHub Actions）

### 6.3 操作步骤

**第一步：创建 GitHub Personal Access Token**

1. 进入 GitHub → Settings → Developer settings → Personal access tokens → Tokens (classic)
2. 点击「Generate new token」
3. 勾选 `repo` 权限（完全控制私有仓库）
4. 生成的 token 请妥善保存

**第二步：使用 API 操作 images.txt**

常用 API 端点：

```bash
# 1. 获取当前 images.txt 内容
curl -H "Authorization: token YOUR_GITHUB_TOKEN" \
     -H "Accept: application/vnd.github.v3+json" \
     https://api.github.com/repos/OWNER/REPO/contents/images.txt

# 2. 获取文件的 SHA（更新时需要）
# 响应中会包含 "sha" 字段，用于更新操作

# 3. 更新 images.txt（添加新镜像）
curl -X PUT -H "Authorization: token YOUR_GITHUB_TOKEN" \
     -H "Accept: application/vnd.github.v3+json" \
     -d '{
       "message": "Add new image",
       "content": "新内容的base64编码",
       "sha": "文件的SHA值"
     }' \
     https://api.github.com/repos/OWNER/REPO/contents/images.txt
```

**示例：Python 脚本添加镜像**

```python
import requests
import base64
import os

GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN")
OWNER = "your-username"
REPO = "docker-mirror"
FILE_PATH = "images.txt"

def get_file_sha():
    url = f"https://api.github.com/repos/{OWNER}/{REPO}/contents/{FILE_PATH}"
    headers = {
        "Authorization": f"token {GITHUB_TOKEN}",
        "Accept": "application/vnd.github.v3+json"
    }
    response = requests.get(url, headers=headers)
    return response.json()["sha"], response.json()["content"]

def add_image(new_image):
    # 获取当前内容和 SHA
    sha, content = get_file_sha()
    current_images = base64.b64decode(content).decode("utf-8")

    # 添加新镜像（如果不存在）
    if new_image not in current_images:
        updated_content = current_images.strip() + f"\n{new_image}\n"

        # 提交更新
        url = f"https://api.github.com/repos/{OWNER}/{REPO}/contents/{FILE_PATH}"
        headers = {
            "Authorization": f"token {GITHUB_TOKEN}",
            "Accept": "application/vnd.github.v3+json"
        }
        data = {
            "message": f"Add image: {new_image}",
            "content": base64.b64encode(updated_content.encode()).decode(),
            "sha": sha
        }
        response = requests.put(url, headers=headers, json=data)
        print(f"Image added: {new_image}")
    else:
        print(f"Image already exists: {new_image}")

# 使用
add_image("ghcr.io/prometheus/prometheus:latest")
```

### 6.4 进阶：定时自动同步

可以结合 GitHub Actions 的定时触发或外部定时任务（如 cron job）来实现自动化：

```yaml
# .github/workflows/auto-sync-list.yml
name: Auto Update Image List
on:
  schedule:
    - cron: '0 0 * * *'  # 每天凌晨执行
  workflow_dispatch:

jobs:
  update-list:
    runs-on: ubuntu-latest
    steps:
      - name: Fetch latest images
        run: |
          # 可以调用其他 API 获取最新镜像版本
          # 然后调用 GitHub API 更新 images.txt
          echo "Fetching latest image versions..."
```
