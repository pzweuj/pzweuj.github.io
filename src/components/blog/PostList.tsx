import Link from 'next/link'
import { BlogPost } from '@/lib/markdown'

export function PostCard({ post }: { post: BlogPost }) {
  return (
    <article className="py-8 space-y-2 group">
      <div className="space-y-3">
        <h2 className="text-2xl font-bold tracking-tight">
          <Link
            href={`/posts/${post.slug}`}
            className="text-gray-900 dark:text-gray-100 hover:text-blue-600 dark:hover:text-blue-400 transition-colors"
          >
            {post.title}
          </Link>
        </h2>
        <div className="flex items-center space-x-4 text-sm text-gray-500 dark:text-gray-400">
          <time dateTime={post.date} className="flex items-center">
            <svg className="mr-1.5 h-4 w-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
              <path 
                strokeLinecap="round" 
                strokeLinejoin="round" 
                strokeWidth={2} 
                d="M8 7V3m8 4V3m-9 8h10M5 21h14a2 2 0 002-2V7a2 2 0 00-2-2H5a2 2 0 00-2 2v12a2 2 0 002 2z" 
              />
            </svg>
            {post.date}
          </time>
        </div>
        <div className="prose dark:prose-invert max-w-none text-gray-600 dark:text-gray-400">
          {post.excerpt}
        </div>
        <div className="flex flex-wrap gap-2">
          {post.tags.map(tag => (
            <Link
              key={tag}
              href={`/tags/${tag}`}
              className="inline-flex items-center px-3 py-1 text-sm bg-gray-100 dark:bg-gray-800 text-gray-700 dark:text-gray-300 rounded-full hover:bg-gray-200 dark:hover:bg-gray-700 transition-colors"
            >
              #{tag}
            </Link>
          ))}
        </div>
      </div>
    </article>
  )
}

export function Pagination({
  currentPage,
  totalPages,
  hasNextPage,
  hasPrevPage,
}: {
  currentPage: number
  totalPages: number
  hasNextPage: boolean
  hasPrevPage: boolean
}) {
  // 生成页码数组
  const getPageNumbers = () => {
    const pages = []
    const showPages = 5 // 显示的页码数量
    
    let start = Math.max(1, currentPage - 2)
    const end = Math.min(totalPages, start + showPages - 1)
    
    if (end - start + 1 < showPages) {
      start = Math.max(1, end - showPages + 1)
    }

    for (let i = start; i <= end; i++) {
      pages.push(i)
    }
    return pages
  }

  return (
    <nav className="flex items-center justify-between border-t border-gray-200 dark:border-gray-800 px-4 sm:px-0 mt-6 py-3">
      <div className="flex flex-1 justify-center gap-1 sm:gap-2">
        {hasPrevPage && (
          <Link
            href={currentPage === 2 ? '/' : `/page/${currentPage - 1}`}
            className="relative inline-flex items-center px-2 sm:px-3 py-1 sm:py-2 text-sm font-medium text-gray-700 dark:text-gray-300 bg-white dark:bg-gray-800 border border-gray-300 dark:border-gray-700 rounded-md hover:bg-gray-50 dark:hover:bg-gray-700 transition-colors"
            aria-label="上一页"
          >
            <svg className="h-5 w-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
            </svg>
          </Link>
        )}

        {getPageNumbers().map(pageNum => (
          <Link
            key={pageNum}
            href={pageNum === 1 ? '/' : `/page/${pageNum}`}
            className={`relative inline-flex items-center px-2 sm:px-4 py-1 sm:py-2 text-sm font-medium rounded-md ${
              pageNum === currentPage
                ? 'bg-blue-500 text-white border border-blue-500'
                : 'text-gray-700 dark:text-gray-300 bg-white dark:bg-gray-800 border border-gray-300 dark:border-gray-700 hover:bg-gray-50 dark:hover:bg-gray-700'
            } transition-colors`}
          >
            {pageNum}
          </Link>
        ))}

        {hasNextPage && (
          <Link
            href={`/page/${currentPage + 1}`}
            className="relative inline-flex items-center px-2 sm:px-3 py-1 sm:py-2 text-sm font-medium text-gray-700 dark:text-gray-300 bg-white dark:bg-gray-800 border border-gray-300 dark:border-gray-700 rounded-md hover:bg-gray-50 dark:hover:bg-gray-700 transition-colors"
            aria-label="下一页"
          >
            <svg className="h-5 w-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
            </svg>
          </Link>
        )}

        <span className="relative inline-flex items-center px-2 sm:px-4 py-1 sm:py-2 text-sm font-medium text-gray-700 dark:text-gray-300">
          {currentPage} / {totalPages} 页
        </span>
      </div>
    </nav>
  )
}