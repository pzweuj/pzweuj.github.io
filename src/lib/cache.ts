import crypto from 'crypto'
import fs from 'fs'
import path from 'path'

const CACHE_DIR = path.join(process.cwd(), '.next', 'cache')
const CACHE_FILE = path.join(CACHE_DIR, 'md-html-cache.json')

export interface CacheEntry {
  hash: string
  html: string
}

export type HtmlCache = Record<string, CacheEntry>

function computeContentHash(pipelineId: string, content: string): string {
  return crypto
    .createHash('md5')
    .update(pipelineId)
    .update(content)
    .digest('hex')
}

export function loadHtmlCache(): HtmlCache {
  try {
    if (fs.existsSync(CACHE_FILE)) {
      return JSON.parse(fs.readFileSync(CACHE_FILE, 'utf-8'))
    }
  } catch {
    // corrupted cache = cache miss
  }
  return {}
}

export function saveHtmlCache(cache: HtmlCache): void {
  try {
    if (!fs.existsSync(CACHE_DIR)) {
      fs.mkdirSync(CACHE_DIR, { recursive: true })
    }
    fs.writeFileSync(CACHE_FILE, JSON.stringify(cache))
  } catch {
    // non-fatal: cache write failure just means cold miss next time
  }
}

export async function renderMarkdownWithCache(
  rawMarkdown: string,
  filePath: string,
  pipelineId: string,
  cache: HtmlCache,
  renderFn: (content: string) => Promise<string>
): Promise<string> {
  const hash = computeContentHash(pipelineId, rawMarkdown)
  const cached = cache[filePath]

  if (cached && cached.hash === hash) {
    return cached.html
  }

  const html = await renderFn(rawMarkdown)
  cache[filePath] = { hash, html }
  return html
}
