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
  const sourcesRef = useRef<Map<HTMLElement, string>>(new Map())
  const renderCountRef = useRef(0)

  // html 变化时重新收集 mermaid 源码（此时占位元素仍在 DOM 中）
  useEffect(() => {
    const container = ref.current
    if (!container) return
    sourcesRef.current.clear()
    container.querySelectorAll('.mermaid').forEach((el) => {
      const source = (el.textContent ?? '').trim()
      if (source) sourcesRef.current.set(el as HTMLElement, source)
    })
  }, [html])

  // 挂载及主题变化时渲染（主题切换时用存储的源码重渲染）
  useEffect(() => {
    renderMermaid()
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [resolvedTheme])

  async function renderMermaid() {
    if (sourcesRef.current.size === 0) return
    const mermaid = (await import('mermaid')).default
    const theme = resolvedTheme === 'dark' ? 'dark' : 'default'
    mermaid.initialize({ startOnLoad: false, securityLevel: 'strict', theme })

    sourcesRef.current.forEach((source, el) => {
      const id = `mermaid-${renderCountRef.current++}`
      mermaid
        .render(id, source)
        .then(({ svg }) => {
          el.innerHTML = svg
        })
        .catch((err) => {
          console.error('Mermaid 渲染失败:', err)
          el.innerHTML = `<pre style="color:#dc2626;white-space:pre-wrap;overflow:auto;">${escapeHtml(source)}</pre>`
        })
    })
  }

  return <div ref={ref} className={className} dangerouslySetInnerHTML={{ __html: html }} />
}