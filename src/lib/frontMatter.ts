import { parse } from 'yaml'

export function parseFrontMatter(source: string, filePath: string) {
  const text = source.replace(/^\uFEFF/, '')
  const firstLine = text.split('\n', 1)[0].trimEnd()
  if (/^---\s*[a-z]/i.test(firstLine) && firstLine !== '---yaml') {
    throw new Error(`Only YAML front matter is supported in ${filePath}`)
  }
  if (firstLine !== '---' && firstLine !== '---yaml') {
    return { data: {}, content: text }
  }

  const rest = text.slice(firstLine.length).replace(/^[ \t]*\r?\n/, '')
  const end = /^---[ \t]*\r?(?:\n|$)/m.exec(rest)
  if (!end) throw new Error(`Unclosed front matter in ${filePath}`)

  const metadata = parse(rest.slice(0, end.index)) ?? {}
  if (typeof metadata !== 'object' || Array.isArray(metadata)) {
    throw new Error(`Expected a YAML mapping in ${filePath}`)
  }
  return {
    data: {
      title: typeof metadata.title === 'string' ? metadata.title : undefined,
      tags: metadata.tags as unknown,
    },
    content: rest.slice(end.index + end[0].length),
  }
}
