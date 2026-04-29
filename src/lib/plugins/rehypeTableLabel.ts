import { visit } from 'unist-util-visit'
import type { Element, Text, Parent } from 'hast'

/**
 * rehype 插件：为表格单元格添加 data-label 属性
 * 用于移动端响应式表格布局
 */
export function rehypeTableLabel() {
  return (tree: Parent) => {
    visit(tree, 'element', (node: Element) => {
      if (node.tagName !== 'table') return

      // 找到所有表头文本
      const headers: string[] = []
      visit(node, 'element', (child: Element) => {
        if (child.tagName === 'th') {
          let text = ''
          visit(child, 'text', (textNode: Text) => {
            text += textNode.value
          })
          headers.push(text.trim())
        }
      })

      if (headers.length === 0) return

      // 为每行的 td 添加 data-label
      visit(node, 'element', (row: Element) => {
        if (row.tagName !== 'tr') return

        const cells = (row.children || []).filter(
          (c) => c.type === 'element' && (c as Element).tagName === 'td'
        ) as Element[]

        cells.forEach((cell: Element, idx: number) => {
          if (idx < headers.length) {
            cell.properties = cell.properties || {}
            cell.properties['data-label'] = headers[idx]
          }
        })
      })
    })
  }
}