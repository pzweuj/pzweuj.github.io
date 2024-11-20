'use client'

import { useState, useMemo } from 'react'
import Link from 'next/link'
import { BlogPost } from '@/lib/markdown'

interface ArchivesListProps {
  initialPosts: BlogPost[]
}

export function ArchivesList({ initialPosts }: ArchivesListProps) {
  const [selectedTag, setSelectedTag] = useState<string | null>(null)

  // 获取所有标签
  const tags = useMemo(() => {
    const tagSet = new Set<string>()
    initialPosts.forEach(post => {
      post.tags.forEach(tag => tagSet.add(tag))
    })
    return Array.from(tagSet).sort()
  }, [initialPosts])

  // 按年份组织文章
  const postsByYear = useMemo(() => {
    // 1. 过滤标签
    const filtered = selectedTag
      ? initialPosts.filter(post => post.tags.includes(selectedTag))
      : initialPosts

    // 2. 先对所有文章按日期从新到旧排序
    const sortedPosts = [...filtered].sort((a, b) => {
      return new Date(b.date).getTime() - new Date(a.date).getTime()
    })

    // 3. 按年份分组
    const grouped: { [key: string]: BlogPost[] } = {}
    sortedPosts.forEach(post => {
      const year = post.date.substring(0, 4)
      if (!grouped[year]) {
        grouped[year] = []
      }
      grouped[year].push(post)
    })

    // 4. 返回按年份排序的数组
    // const years = Object.keys(grouped).sort((a, b) => Number(a) - Number(b))
    const years = Object.keys(grouped).sort((a, b) => Number(b) - Number(a))

    // 创建年份和文章的有序数组
    const sortedEntries = years.map(year => ({
      year,
      posts: grouped[year]
    }))

    return sortedEntries
  }, [initialPosts, selectedTag])

  return (
    <div>
      {/* 标签过滤器 */}
      <div className="mb-8">
        <div className="flex flex-wrap gap-2">
          <button
            onClick={() => setSelectedTag(null)}
            className={`px-3 py-1 rounded-full text-sm transition-colors ${
              selectedTag === null
                ? 'bg-blue-500 text-white'
                : 'bg-gray-100 dark:bg-gray-800 text-gray-700 dark:text-gray-300 hover:bg-gray-200 dark:hover:bg-gray-700'
            }`}
          >
            全部
          </button>
          {tags.map(tag => (
            <button
              key={tag}
              onClick={() => setSelectedTag(tag)}
              className={`px-3 py-1 rounded-full text-sm transition-colors ${
                selectedTag === tag
                  ? 'bg-blue-500 text-white'
                  : 'bg-gray-100 dark:bg-gray-800 text-gray-700 dark:text-gray-300 hover:bg-gray-200 dark:hover:bg-gray-700'
              }`}
            >
              {tag}
            </button>
          ))}
        </div>
      </div>

      {/* 文章列表 */}
      <div className="space-y-12">
        {postsByYear.map(({ year, posts }) => (
          <div key={year}>
            <h2 className="text-2xl font-bold mb-4 text-gray-900 dark:text-gray-100">
              {year}
              <span className="ml-2 text-sm text-gray-500">
                ({posts.length} 篇)
              </span>
            </h2>
            <div className="space-y-4">
              {posts.map(post => (
                <article
                  key={post.slug}
                  className="flex items-baseline gap-4 group"
                >
                  <time className="w-24 text-sm text-gray-500 dark:text-gray-400">
                    {post.date.substring(5)}
                  </time>
                  <h3 className="flex-1">
                    <Link
                      href={`/posts/${post.slug}`}
                      className="text-gray-900 dark:text-gray-100 hover:text-blue-600 dark:hover:text-blue-400 transition-colors"
                    >
                      {post.title}
                    </Link>
                  </h3>
                  <div className="hidden md:flex gap-2">
                    {post.tags.map(tag => (
                      <button
                        key={tag}
                        onClick={() => setSelectedTag(tag)}
                        className="text-sm text-gray-500 dark:text-gray-400 hover:text-blue-600 dark:hover:text-blue-400"
                      >
                        #{tag}
                      </button>
                    ))}
                  </div>
                </article>
              ))}
            </div>
          </div>
        ))}
      </div>

      {/* 无结果提示 */}
      {postsByYear.length === 0 && (
        <div className="text-center py-12 text-gray-500 dark:text-gray-400">
          没有找到相关文章
        </div>
      )}
    </div>
  )
} 