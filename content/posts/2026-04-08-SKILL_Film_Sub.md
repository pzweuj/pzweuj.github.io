---
title: 电影记录和订阅管理
tags: software
---

> 给自己的小龙虾和爱马仕快速开发了两个服务



## Sub Manager Lite

订阅管理服务，通过自然语言管理你的订阅。

### 部署

**Docker（推荐）**

```yaml
services:
  sub-manager:
    image: ghcr.io/pzweuj/sub_manager_lite:latest
    ports:
      - "8000:8000"
    environment:
      - API_TOKEN=your-secret-token
      - BASE_CURRENCY=CNY
    volumes:
      - ./data:/data
```

```bash
docker-compose up -d
```

**本地运行**

```bash
pip install -r requirements.txt
export API_TOKEN="your-token"
uvicorn app.main:app --reload
```

### 使用

将 [SKILL.md](https://github.com/pzweuj/sub_manager_lite/raw/refs/heads/main/SKILL.md) 地址给你的 AI Agent（OpenClaw、Hermes-Agent 等），首次使用时提供服务 URL 和 Token 即可。

**示例对话**

```
用户："帮我记录一个订阅：ChatGPT Plus，每月 20 美元，到期日 2025-01-15"
Agent：✓ 已创建订阅

用户："今年订阅总共花了多少？"
Agent：年度总花费 ¥1,650.00（生产力 ¥1,440 + 娱乐 ¥210）

用户："最近有什么订阅要扣费？"
Agent：未来 7 天：Netflix $15.99、ChatGPT Plus $20.00、GitHub Copilot $10.00

用户："把我的 Wallos 数据迁移过来"
Agent：请提供 Wallos 的 URL 和 API Token
用户：URL 是 http://wallos.example.com，Token 是 xxx
Agent：✓ 已从 Wallos 迁移 12 个订阅记录
```

**GitHub**: [https://github.com/pzweuj/sub_manager_lite](https://github.com/pzweuj/sub_manager_lite)



## Film Record Lite

观影记录服务，通过自然语言管理你看过的电影。

### 部署

**Docker（推荐）**

```yaml
services:
  film-record-lite:
    image: ghcr.io/pzweuj/film_record_lite:latest
    ports:
      - "8000:8000"
    environment:
      - FILM_RECORD_TOKEN=your-secret-token
    volumes:
      - ./data:/app/data
```

```bash
docker compose up -d
```

**本地运行**

```bash
pip install -e .
FILM_RECORD_TOKEN=your_token python -m film_record_lite.server --port 8000
```

### 使用

将 [SKILL.md](https://github.com/pzweuj/film_record_lite/raw/refs/heads/main/SKILL.md) 地址给你的 AI Agent，首次使用时提供服务 URL 和 Token 即可。

**示例对话**

```
用户："我看了无间道，评分 9 分"
Agent：✓ 已记录「无间道」主演：刘德华,梁朝伟,黄秋生,曾志伟，评分：9/10

用户："我最近看了无间道 9 分、肖申克的救赎满分、霸王别姬 9.5 分"
Agent：✓ 已添加 3 部电影记录

用户："我看过哪些梁朝伟的电影？"
Agent：「梁朝伟」主演的电影：无间道（评分 10/10）

用户："把无间道的评分改成 10 分"
Agent：✓ 已更新评分
```

**GitHub**: [https://github.com/pzweuj/film_record_lite](https://github.com/pzweuj/film_record_lite)