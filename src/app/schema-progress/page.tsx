import { getAllSchemaProgress } from '@/lib/markdown'
import Link from 'next/link'

export default async function SchemaProgressPage() {
  const posts = await getAllSchemaProgress()

  // 按年份组织文章
  const postsByYear: { [key: string]: typeof posts } = {}
  posts.forEach(post => {
    const year = post.date.substring(0, 4)
    if (!postsByYear[year]) {
      postsByYear[year] = []
    }
    postsByYear[year].push(post)
  })

  const years = Object.keys(postsByYear).sort((a, b) => Number(b) - Number(a))

  return (
    <div className="max-w-4xl mx-auto px-4 py-12">
      <h1 className="text-3xl font-bold mb-8 text-gray-900 dark:text-gray-100">
        Schema进度
      </h1>

      <p className="text-gray-600 dark:text-gray-400 mb-8">
        记录 <a href="https://github.com/SchemaBio/schema-platform" className="text-blue-600 dark:text-blue-400 hover:underline" target="_blank" rel="noopener noreferrer">schema-platform</a> 项目开发进度。
      </p>

      {years.length === 0 ? (
        <div className="text-center py-12 text-gray-500 dark:text-gray-400">
          暂无进度记录
        </div>
      ) : (
        <div className="space-y-12">
          {years.map(year => (
            <div key={year}>
              <h2 className="text-2xl font-bold mb-4 text-gray-900 dark:text-gray-100">
                {year}
                <span className="ml-2 text-sm text-gray-500">
                  ({postsByYear[year].length} 篇)
                </span>
              </h2>
              <div className="space-y-4">
                {postsByYear[year].map(post => (
                  <article
                    key={post.slug}
                    className="flex items-baseline gap-4 group"
                  >
                    <time className="w-24 text-sm text-gray-500 dark:text-gray-400">
                      {post.date.substring(5)}
                    </time>
                    <h3 className="flex-1">
                      <Link
                        href={`/schema-progress/${post.slug}`}
                        className="text-gray-900 dark:text-gray-100 hover:text-blue-600 dark:hover:text-blue-400 transition-colors"
                      >
                        {post.title}
                      </Link>
                    </h3>
                  </article>
                ))}
              </div>
            </div>
          ))}
        </div>
      )}
    </div>
  )
}