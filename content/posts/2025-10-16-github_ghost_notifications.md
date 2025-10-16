---
title: 删除Github的幽灵通知
tags: coding
---

最近github的通知里一直都有小蓝点提示有未读的通知，但是点进去后又啥都没有（但是可以看到通知来源的仓库）。

查了一下，问题是这些仓库是钓鱼的，他们通过在仓库里@，然后github的官方就会用官方的noreply邮箱给你发邮件，让你误认为是真的。

github官方把这些仓库删了，但留下了小蓝点通知的手尾。

### 下面是删除这些幽灵通知的方案

1. 点击自己头像 -> Setting -> Developer Settings -> Personal access tokens -> Tokens (classic) -> Gernerate new token ->  Gernerate new token(classic) 创建一个包含了 `notifications`权限的token。

**这个token如果以后不再用，就设定最短的期限，或者用完即删！**

2. 在linux中使用下面的脚本即可完成删除幽灵通知

```bash
#!/bin/bash

# 检查jq是否安装
if ! command -v jq &> /dev/null
then
    echo "错误：jq 未安装。请先安装 jq (例如：sudo apt-get install jq 或 brew install jq)。"
    exit 1
fi

echo "GitHub 未读通知清理工具 (修复 204 成功检测)"
echo "--------------------------"

# 安全地获取GitHub个人访问令牌
read -s -p "请输入你的GitHub个人访问令牌 (ghp_XXX)：" GITHUB_TOKEN
echo "" # 换行

if [ -z "$GITHUB_TOKEN" ]; then
    echo "错误：未输入令牌。无法继续。"
    exit 1
fi

BASE_URL="https://api.github.com"
HEADERS=(
  -H "Accept: application/vnd.github+json"
  -H "Authorization: token $GITHUB_TOKEN"
  -H "X-GitHub-Api-Version: 2022-11-28"
)

echo "正在获取未读通知..."

# 获取所有未读通知
NOTIFICATIONS_GET_RESPONSE=$(curl -L -s -w "%{http_code}" "${HEADERS[@]}" "$BASE_URL/notifications?participating=true")

HTTP_STATUS_GET=$(echo "$NOTIFICATIONS_GET_RESPONSE" | tail -n1)
NOTIFICATIONS_JSON=$(echo "$NOTIFICATIONS_GET_RESPONSE" | sed '$d') # 移除状态码以获取实际响应体

echo "GET /notifications HTTP状态码: $HTTP_STATUS_GET"

# 检查获取通知的HTTP状态码
if [ "$HTTP_STATUS_GET" -ne 200 ] && [ "$HTTP_STATUS_GET" -ne 204 ]; then # 200 OK, 204 No Content (No Notificatons)
    echo "错误：获取通知失败，HTTP状态码 $HTTP_STATUS_GET。"
    echo "原始响应：$NOTIFICATIONS_JSON"
    if echo "$NOTIFICATIONS_JSON" | grep -q "Bad credentials"; then
        echo "API返回 'Bad credentials'。你的个人访问令牌可能不正确或已过期/撤销，或者权限不足。"
    fi
    exit 1
fi

# 如果返回204，表示没有内容，即没有未读通知
if [ "$HTTP_STATUS_GET" -eq 204 ]; then
    echo "API返回 204 (No Content)。这意味着目前没有未读通知。"
    exit 0
fi

# 使用jq解析通知ID
NOTIFICATION_IDS=$(echo "$NOTIFICATIONS_JSON" | jq -r '.[].id' 2>/dev/null) # 2>/dev/null 隐藏jq错误信息

# 检查jq是否成功提取ID
if [ $? -ne 0 ]; then
    echo "错误：jq解析JSON失败。请确认API响应是有效的JSON格式。"
    echo "原始响应：$NOTIFICATIONS_JSON"
    exit 1
fi

if [ -z "$NOTIFICATION_IDS" ]; then
    echo "没有找到可提取的未读通知ID。"
    echo "----------------------------------------------------"
    echo "以下是GET /notifications的原始JSON响应（供调试）："
    echo "----------------------------------------------------"
    echo "$NOTIFICATIONS_JSON"
    echo "----------------------------------------------------"
    echo "请检查以上JSON，看是否有预期的通知对象和它们各自的'id'字段。"
else
    echo "找到以下未读通知ID："
    echo "$NOTIFICATION_IDS"
    echo "即将逐个删除这些通知..."

    for ID in $NOTIFICATION_IDS; do
        DELETE_URL="$BASE_URL/notifications/threads/$ID"
        echo "  - 正在尝试删除通知 ID: $ID ..."
        
        DELETE_RESPONSE=$(curl -L -s -X DELETE "${HEADERS[@]}" "$DELETE_URL" -w "%{http_code}")
        HTTP_STATUS_DELETE=$(echo "$DELETE_RESPONSE" | tail -n1)
        RESPONSE_BODY=$(echo "$DELETE_RESPONSE" | sed '$d') # 移除状态码以获取实际响应体

        # === 关键修改：将 204 也视为成功 ===
        if [ "$HTTP_STATUS_DELETE" -eq 205 ] || [ "$HTTP_STATUS_DELETE" -eq 204 ]; then
            echo "    ✅ 通知 ID: $ID 已成功标记为已读 (HTTP $HTTP_STATUS_DELETE)."
        elif [ "$HTTP_STATUS_DELETE" -eq 404 ]; then
            echo "    ⚠️ 通知 ID: $ID 未找到或已被处理 (HTTP 404)。"
        elif [ "$HTTP_STATUS_DELETE" -eq 401 ]; then
            echo "    ❌ 错误：删除通知 ID: $ID 时遇到 'Bad credentials' (HTTP 401)。请检查令牌权限。"
            echo "       响应：$RESPONSE_BODY"
            break # 令牌有问题，停止后续删除
        else
            echo "    ❌ 错误：删除通知 ID: $ID 失败 (HTTP $HTTP_STATUS_DELETE)。"
            echo "       响应：$RESPONSE_BODY"
        fi
    done
    echo "--------------------------"
    echo "所有通知处理完毕。"
fi

exit 0
```

