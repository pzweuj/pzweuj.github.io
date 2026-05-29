import { visit } from 'unist-util-visit'
import type { Element, Parent } from 'hast'

/**
 * rehype 插件：为代码块添加一键复制按钮
 * 在 pre > code 结构的 pre 外包裹 div.code-block-wrapper
 * 仅做 AST 变换，事件委托由 CodeCopyScript 组件负责
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

      const button: Element = {
        type: 'element',
        tagName: 'button',
        properties: {
          className: ['code-copy-btn'],
          'aria-label': '复制代码',
          'data-action': 'copy-code',
        },
        children: [],
      }

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
