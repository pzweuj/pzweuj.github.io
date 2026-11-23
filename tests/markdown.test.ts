import assert from 'node:assert/strict'
import test from 'node:test'
import { createMarkdownProcessor } from '../src/lib/createMarkdownProcessor'
import { parseFrontMatter } from '../src/lib/frontMatter'

test('reads YAML and preserves posts without metadata', () => {
  const body = '# 正文\r\n'
  for (const header of ['', '---\r\ntitle: NIPT\r\ntags: [coding]\r\n---\r\n']) {
    const post = parseFrontMatter(header + body, 'post.md')
    assert.equal(post.content, body)
    if (header) assert.deepEqual(post.data, { title: 'NIPT', tags: ['coding'] })
  }
  assert.throws(() => parseFrontMatter('---js\ntitle: sample\n---\n', 'post.md'), /Only YAML/)
})

test('filters author-supplied frames and preserves supported Markdown features', async () => {
  const markdown = [
    '<iframe src="https://example.org/"></iframe>',
    '[qqmusic:123]',
    '<font color="red">格式</font>',
    '$x^2$',
    '```mermaid\ngraph TD\nA --> B\n```',
    '```js\nconst x = 1\n```',
    '[下载](ftp://example.org/file)',
    '注释[^1]',
    '[^1]: 脚注',
  ].join('\n\n')
  const processor = createMarkdownProcessor({ enableCodeCopy: true })
  const html = String(await processor.process(markdown))
  assert.equal(html.match(/<iframe\b/g)?.length, 1)
  assert.match(html, /src="https:\/\/i\.y\.qq\.com\/[^"]+songid=123"/)
  assert.match(html, /<font color="red">格式<\/font>/)
  assert.match(html, /class="katex"/)
  assert.match(html, /class="mermaid"/)
  assert.match(html, /class="code-copy-btn"/)
  assert.match(html, /href="ftp:\/\/example\.org\/file"/)
  assert.match(html, /href="#user-content-fn-1"/)
  assert.match(html, /id="user-content-fn-1"/)
})
