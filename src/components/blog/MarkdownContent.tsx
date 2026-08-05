'use client'

import { useEffect, useRef } from 'react'
import { useTheme } from 'next-themes'

function escapeHtml(s: string) {
  return s
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
}

/**
 * 在浏览器端渲染 Markdown 内容，并负责把 mermaid 占位元素渲染成 SVG。
 * 服务端只输出 <div class="mermaid">源码</div>，实际的图形渲染依赖浏览器，
 * 因此这里必须是客户端组件。
 *
 * mermaid 体积较大，这里用动态 import 按需加载：只有页面里确实存在
 * .mermaid 块时才拉取，普通文章不会白白加载这份 JS。
 */
export default function MarkdownContent({
  html,
  className = 'prose max-w-none',
}: {
  html: string
  className?: string
}) {
  const ref = useRef<HTMLDivElement>(null)
  const { resolvedTheme } = useTheme()
  // 渲染序号：每次渲染递增，异步完成时若发现已有更新的渲染，则丢弃本次结果。
  // 主题切换时（尤其首次挂载 resolvedTheme 为 undefined 后解析成真实主题）会
  // 并发触发多次 renderMermaid，不丢弃旧结果的话，慢的旧主题可能覆盖新主题。
  const renderSeqRef = useRef(0)
  // mermaid 需要唯一的 diagram id，重复 id 会抛错，用全局递增计数器保证唯一。
  const idCounterRef = useRef(0)

  // 只在主题确定后渲染；主题或 html 变化时重渲染（从 data 属性重新取源码）。
  useEffect(() => {
    if (!resolvedTheme) return
    const container = ref.current
    if (!container) return
    renderMermaid(container)
  }, [resolvedTheme, html])

  async function renderMermaid(container: HTMLElement) {
    const blocks = Array.from(container.querySelectorAll<HTMLElement>('.mermaid'))
    if (blocks.length === 0) return

    const mermaid = (await import('mermaid')).default
    const theme = resolvedTheme === 'dark' ? 'dark' : 'default'
    mermaid.initialize({ startOnLoad: false, securityLevel: 'strict', theme })

    const seq = ++renderSeqRef.current
    for (const el of blocks) {
      // 源码优先从 data-mermaid-src 读取（首次渲染后 textContent 已被 SVG 覆盖）
      const source = el.getAttribute('data-mermaid-src') ?? (el.textContent ?? '').trim()
      if (!source) continue

      const id = `mermaid-${idCounterRef.current++}`
      try {
        const { svg } = await mermaid.render(id, source)
        if (seq !== renderSeqRef.current) return // 已有更新的渲染，丢弃过期结果
        el.innerHTML = svg
      } catch (err) {
        console.error('Mermaid 渲染失败:', err)
        if (seq === renderSeqRef.current) {
          el.innerHTML = `<pre style="color:#dc2626;white-space:pre-wrap;overflow:auto;">${escapeHtml(source)}</pre>`
        }
      }
    }
  }

  return <div ref={ref} className={className} dangerouslySetInnerHTML={{ __html: html }} />
}