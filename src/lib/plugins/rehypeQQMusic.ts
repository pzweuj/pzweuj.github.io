import { visit } from 'unist-util-visit'
import type { Root } from 'hast'

// Only this shortcode can create an iframe, after author HTML has been sanitized.
export function rehypeQQMusic() {
  return (tree: Root) => {
    visit(tree, 'element', node => {
      if (node.tagName !== 'p' || node.children.length !== 1) return
      const child = node.children[0]
      if (child.type !== 'text') return
      const song = /^\[qqmusic:(\d+)\]$/.exec(child.value)
      if (!song) return

      node.tagName = 'div'
      node.properties = { className: ['qq-music-player'] }
      node.children = [{
        type: 'element',
        tagName: 'iframe',
        properties: {
          title: 'QQ音乐',
          src: `https://i.y.qq.com/n2/m/outchain/player/index.html?songid=${song[1]}`,
          width: '100%',
          height: 86,
          frameBorder: 0,
        },
        children: [],
      }]
    })
  }
}
