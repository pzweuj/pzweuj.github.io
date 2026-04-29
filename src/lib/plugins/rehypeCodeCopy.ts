import { visit } from 'unist-util-visit'
import type { Element, Parent } from 'hast'

/**
 * rehype 插件：为代码块添加一键复制按钮
 * 在 pre > code 结构的 pre 外包裹 div.code-block-wrapper
 */
export function rehypeCodeCopy() {
  return (tree: Parent) => {
    visit(tree, 'element', (node: Element, index, parent) => {
      if (node.tagName !== 'pre') return
      if (!parent || typeof index !== 'number') return

      const codeEl = (node.children || []).find(
        (c) => c.type === 'element' && (c as Element).tagName === 'code'
      ) as Element | undefined
      if (!codeEl) return

      // 复制按钮（VitePress 风格：纯图标 + 复制后左侧滑出"已复制"标签）
      const button: Element = {
        type: 'element',
        tagName: 'button',
        properties: {
          className: ['code-copy-btn'],
          'aria-label': '复制代码',
          onclick: "(function(btn){" +
            "var code=btn.parentElement.querySelector('code');" +
            "if(!code)return;" +
            "var t=code.innerText;" +
            "var done=function(){" +
              "btn.classList.add('copied');" +
              "setTimeout(function(){btn.classList.remove('copied')},2000);" +
            "};" +
            "if(navigator.clipboard){navigator.clipboard.writeText(t).then(done)}" +
            "else{var ta=document.createElement('textarea');ta.value=t;document.body.appendChild(ta);ta.select();document.execCommand('copy');document.body.removeChild(ta);done()}" +
          "})(this)",
        },
        children: [],
      }

      // 包裹层
      const wrapper: Element = {
        type: 'element',
        tagName: 'div',
        properties: { className: ['code-block-wrapper'] },
        children: [button, node],
      }

      const parentEl = parent as Element
      parentEl.children.splice(index, 1, wrapper)
    })
  }
}
