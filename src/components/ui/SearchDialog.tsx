'use client'

import { useEffect, useState } from 'react'
import { useRouter } from 'next/navigation'
import { SearchIcon } from './Icons'
import { Command } from 'cmdk'
import type { SearchResult } from '@/lib/search'

export default function SearchDialog() {
  const [open, setOpen] = useState(false)
  const [query, setQuery] = useState('')
  const [searchIndex, setSearchIndex] = useState<SearchResult[]>([])
  const [results, setResults] = useState<SearchResult[]>([])
  const router = useRouter()

  // 加载搜索索引
  useEffect(() => {
    fetch('/api/search')
      .then(res => res.json())
      .then(data => setSearchIndex(data))
      .catch(error => console.error('Failed to load search index:', error))
  }, [])

  // 本地搜索处理
  useEffect(() => {
    if (!query.trim()) {
      setResults([])
      return
    }

    const searchQuery = query.toLowerCase()
    const filtered = searchIndex.filter(item => 
      item.title.toLowerCase().includes(searchQuery) ||
      item.excerpt.toLowerCase().includes(searchQuery)
    )

    const articleMap = new Map()
    filtered.forEach((item: SearchResult) => {
      const existing = articleMap.get(item.slug)
      if (!existing || (item.matchType === 'title' && existing.matchType === 'content')) {
        articleMap.set(item.slug, item)
      }
    })

    setResults(Array.from(articleMap.values()))
  }, [query, searchIndex])

  // 键盘快捷键处理
  useEffect(() => {
    const down = (e: KeyboardEvent) => {
      if (e.key === 'k' && (e.metaKey || e.ctrlKey)) {
        e.preventDefault()
        setOpen((open) => !open)
      }
    }

    document.addEventListener('keydown', down)
    return () => document.removeEventListener('keydown', down)
  }, [])

  return (
    <>
      <button
        onClick={() => setOpen(true)}
        className="group p-2 hover:bg-gray-100 dark:hover:bg-gray-800 rounded-full transition-colors"
        aria-label="搜索 (⌘K)"
      >
        <SearchIcon className="w-5 h-5 text-gray-600 dark:text-gray-400 group-hover:text-gray-900 dark:group-hover:text-gray-100" />
      </button>

      {open && (
        <div className="fixed inset-0 z-50">
          <div className="fixed inset-0 bg-gray-900/50 backdrop-blur-sm" onClick={() => setOpen(false)} />
          
          <div className="fixed inset-x-4 top-8 mx-auto max-w-2xl">
            <div className="relative bg-white dark:bg-gray-900 rounded-xl shadow-2xl ring-1 ring-gray-200 dark:ring-gray-800">
              <Command className="relative" loop>
                <div className="flex items-center px-4 border-b border-gray-200 dark:border-gray-800">
                  <SearchIcon className="w-5 h-5 text-gray-400" />
                  <Command.Input
                    value={query}
                    onValueChange={setQuery}
                    placeholder="搜索文章..."
                    className="flex-1 h-14 px-4 font-medium bg-transparent outline-none placeholder:text-gray-400"
                  />
                  {query && (
                    <button
                      onClick={() => setQuery('')}
                      className="text-xs text-gray-400 hover:text-gray-600 dark:hover:text-gray-300"
                    >
                      清除
                    </button>
                  )}
                </div>

                <Command.List className="max-h-[60vh] overflow-y-auto overscroll-contain p-2">
                  {query && results.length === 0 ? (
                    <Command.Empty className="py-6 text-center text-sm text-gray-500">
                      未找到相关文章
                    </Command.Empty>
                  ) : (
                    results.map((result) => (
                      <Command.Item
                        key={result.slug}
                        value={result.title}
                        onSelect={() => {
                          router.push(`/posts/${result.slug}`)
                          setOpen(false)
                        }}
                        className="group px-4 py-3 rounded-lg cursor-pointer data-[selected=true]:bg-gray-100 dark:data-[selected=true]:bg-gray-800"
                      >
                        <div className="space-y-1">
                          <h3 className="font-medium text-gray-900 dark:text-gray-100 line-clamp-1">
                            {result.title}
                          </h3>
                          <p className="text-sm text-gray-500 dark:text-gray-400 line-clamp-2">
                            {result.excerpt}
                          </p>
                          <div className="flex items-center gap-2 text-xs text-gray-400">
                            <time>{result.date}</time>
                          </div>
                        </div>
                      </Command.Item>
                    ))
                  )}
                </Command.List>

                {results.length > 0 && query && (
                  <div className="flex items-center h-10 px-4 border-t border-gray-200 dark:border-gray-800">
                    <p className="flex-1 text-xs text-gray-400">
                      <kbd className="px-1.5 py-0.5 text-xs rounded bg-gray-100 dark:bg-gray-800">↵</kbd> 选择
                    </p>
                    <p className="text-xs text-gray-400">找到 {results.length} 篇文章</p>
                  </div>
                )}
              </Command>
            </div>
          </div>
        </div>
      )}
    </>
  )
} 