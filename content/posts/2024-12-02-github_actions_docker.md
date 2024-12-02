---
title: 使用Github Actions进行docker build
tags: coding
---

## 效果

建立了一个仓库，基本功能是当检测到仓库中存在dockerfile更新时，自动进行docker build。

镜像名称会以子文件夹的名称进行命名，而tag则会以git commit该dockerfile更新的comment设定。

另外，稍微改造一下这个方案，就可以实现其他项目的自动docker build了。


示例仓库：[MyDockerImagePublic](https://github.com/pzweuj/MyDockerImagePublic)

构建结果：[Packages](https://github.com/pzweuj?tab=packages)

例如：

```
- repo
    - sub_01
        - dockerfile
    - sub_02
        - dockerfile
- dockerfile
```



在上面这种结构下，会建立

```
ghcr.io/pzweuj:commit1
ghcr.io/pzweuj/sub_01:commit2
ghcr.io/pzweuj/sub_02:commit3
```



## 配置

将下面的yml配置文件放置在仓库的

```
.github/workflows
```

文件夹下

配置文件参考

```yaml
name: Build and Push Docker Images

on:
  push:
    paths:
      - '**/[Dd]ockerfile'
      - '**/Dockerfile*'

jobs:
  detect-changes:
    runs-on: ubuntu-latest
    outputs:
      changed_files: ${{ steps.detect.outputs.changed_files }}
    steps:
      - name: Checkout code
        uses: actions/checkout@v3
        with:
          fetch-depth: 0

      - name: Detect changed Dockerfiles
        id: detect
        run: |
          echo "Current SHA: ${{ github.sha }}"
          echo "Previous SHA: ${{ github.event.before }}"
          
          # 显示所有变更的文件
          echo "All changed files:"
          git diff-tree --no-commit-id --name-only -r ${{ github.sha }}
          
          # 尝试使用不同的git命令来检测变更
          CHANGED_FILES=$(git diff-tree --no-commit-id --name-only -r ${{ github.sha }} | grep -i dockerfile || true)
          echo "changed_files=${CHANGED_FILES}" >> $GITHUB_OUTPUT

      - name: Extract Commit Message
        id: extract_message
        run: |
          COMMIT_MSG=$(git log -1 --pretty=%B)
          FORMATTED_MSG=$(echo "$COMMIT_MSG" | tr ' ' '-' | tr -d '\n' | tr -cd '[:alnum:]-')
          echo "Commit Message: $FORMATTED_MSG"
          echo "COMMIT_MESSAGE=$FORMATTED_MSG" >> $GITHUB_ENV

      - name: Debug Output
        run: |
          echo "GITHUB_OUTPUT content:"
          cat $GITHUB_OUTPUT
          echo "Changed files value:"
          echo "${{ steps.detect.outputs.changed_files }}"

  build:
    needs: detect-changes
    if: ${{ contains(needs.detect-changes.outputs.changed_files, 'Dockerfile') }}
    runs-on: ubuntu-latest
    steps:
      - name: Checkout code
        uses: actions/checkout@v3

      - name: Login to GitHub Container Registry
        uses: docker/login-action@v2
        with:
          registry: ghcr.io
          username: ${{ github.repository_owner }}
          password: ${{ secrets.GHCR_PAT }}

      - name: Build and push images
        run: |
          # 获取 commit message 并清理特殊字符
          VERSION="${{ github.event.head_commit.message }}"
          # 移除特殊字符，只保留字母、数字、点、横杠和下划线
          VERSION=$(echo "$VERSION" | tr -dc 'a-zA-Z0-9._-')
          
          echo "Using version tag: $VERSION"
          
          for file in ${{ needs.detect-changes.outputs.changed_files }}; do
            echo "Building and pushing image for $file..."
            dir=$(dirname $file)
            docker build -t ghcr.io/pzweuj/$dir:$VERSION -t ghcr.io/pzweuj/$dir:latest -f $file .
            echo "Pushing image ghcr.io/pzweuj/$dir:$VERSION"
            docker push ghcr.io/pzweuj/$dir:$VERSION
            docker push ghcr.io/pzweuj/$dir:latest
          done

      - name: Debug Changed Files
        run: |
          echo "Changed files from previous job:"
          echo "${{ needs.detect-changes.outputs.changed_files }}"
```



## 环境变量

需要对这个仓库进行环境变量设置，因为我们需要将镜像推送到ghcr中，因此需要配置token。在仓库的设置中进行配置，避免token泄漏。



### token获取

点击自己的头像 -> Settings -> Developer settings -> Personal access tokens -> Tokens(classic)  -> Generate new token -> Generate new token(classic)

在验证账户后，仅勾选 “读”、“写”权限，token有效期最长可选择一年，生成后只会在页面显示一次，记得复制下来。



### 仓库配置

在仓库页面，点击上方的 Settings，然后选择左侧的 Secrets and variables -> Actions，添加一个Secrets，其中名称是GHCR\_PAT，值是刚刚获取的token。保存后，这个仓库即可实现对传输进的dockerfile进行自动docker build的功能。









