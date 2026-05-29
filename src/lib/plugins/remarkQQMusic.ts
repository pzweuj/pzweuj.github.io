import { visit, SKIP } from 'unist-util-visit'
import type { Plugin } from 'unified'
import type { Root, Paragraph, Text, HTML } from 'mdast'

interface QQMusicOptions {
  className?: string
}

const QQ_MUSIC_REGEX = /^\[qqmusic:(\d+)\]$/

export const remarkQQMusic: Plugin<[QQMusicOptions?], Root> = (options = {}) => {
  const className = options.className || 'qq-music-player'

  return (tree) => {
    visit(tree, 'paragraph', (node: Paragraph, index, parent) => {
      if (
        node.children.length === 1 &&
        node.children[0].type === 'text'
      ) {
        const textNode = node.children[0] as Text
        const match = textNode.value.match(QQ_MUSIC_REGEX)
        if (match && typeof index === 'number' && parent) {
          const songId = match[1]
          const htmlNode: HTML = {
            type: 'html',
            value: `<div class="${className}"><iframe frameborder="no" border="0" marginwidth="0" marginheight="0" width="100%" height="86" src="https://i.y.qq.com/n2/m/outchain/player/index.html?songid=${songId}"></iframe></div>`,
          }
          parent.children.splice(index, 1, htmlNode)
          return [SKIP, index]
        }
      }
    })
  }
}
