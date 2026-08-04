import { visit } from 'unist-util-visit'
import type { Element, Parent, Text } from 'hast'

const MERMAID_CLASS = 'language-mermaid'

/**
 * rehype 插件：把 mermaid 代码块转换为客户端可渲染的占位元素。
 *
 * mermaid 需要浏览器运行时才能把图形渲染成 SVG，因此这里不做服务端渲染，
 * 而是把:
 *   <pre><code class="language-mermaid">…源码…</code></pre>
 * 替换为:
 *   <div class="mermaid">…源码…</div>
 * 由客户端组件 MarkdownContent 在挂载后调用 mermaid.render() 完成渲染。
 *
 * ⚠️ 必须在 rehype-prism 之前运行，否则 prism 会把 mermaid 源码当作普通
 * 代码高亮并加上行号，从而破坏图形源码。
 */
export function rehypeMermaid() {
  return (tree: Parent) => {
    visit(tree, 'element', (node: Element, index, parent) => {
      if (node.tagName !== 'pre') return

      const code = node.children?.find(
        (c) => c.type === 'element' && (c as Element).tagName === 'code'
      ) as Element | undefined
      if (!code) return

      const className = code.properties?.className
      const classes = Array.isArray(className) ? className.map(String) : []
      if (!classes.includes(MERMAID_CLASS)) return

      const grandparent = parent as Parent | undefined
      if (!grandparent) return

      // 收集源码文本（mermaid 源码可能被拆成多个文本节点）
      let source = ''
      const collect = (n: Text | Element) => {
        if (n.type === 'text') source += (n as Text).value
        else if (n.children) (n.children as (Text | Element)[]).forEach(collect)
      }
      ;(code.children as (Text | Element)[]).forEach(collect)
      source = source.trim()

      const mermaidDiv: Element = {
        type: 'element',
        tagName: 'div',
        properties: { className: ['mermaid'] },
        children: [{ type: 'text', value: source }],
      }

      const idx = grandparent.children.indexOf(node)
      if (idx !== -1) {
        grandparent.children[idx] = mermaidDiv
      }
    })
  }
}