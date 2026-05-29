'use client'

import { useEffect } from 'react'

/**
 * 代码块复制按钮的事件委托。
 * 配合 rehypeCodeCopy 插件使用——插件负责 DOM 结构，本组件负责交互。
 */
export function CodeCopyScript() {
  useEffect(() => {
    function handleClick(e: MouseEvent) {
      const btn = (e.target as HTMLElement).closest?.('[data-action="copy-code"]')
      if (!btn) return

      const code = btn.parentElement?.querySelector('code')
      if (!code) return

      const text = code.innerText
      const done = () => {
        btn.classList.add('copied')
        setTimeout(() => btn.classList.remove('copied'), 2000)
      }

      if (navigator.clipboard) {
        navigator.clipboard.writeText(text).then(done)
      } else {
        const ta = document.createElement('textarea')
        ta.value = text
        document.body.appendChild(ta)
        ta.select()
        document.execCommand('copy')
        document.body.removeChild(ta)
        done()
      }
    }

    document.addEventListener('click', handleClick)
    return () => document.removeEventListener('click', handleClick)
  }, [])

  return null
}
