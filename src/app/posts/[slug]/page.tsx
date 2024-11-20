import { getAllPosts } from '@/lib/markdown'
import { notFound } from 'next/navigation'
import { Metadata } from 'next'

interface Props {
  params: { slug: string }
}

// 生成动态元数据
export async function generateMetadata({ params }: Props): Promise<Metadata> {
  const posts = await getAllPosts()
  const post = posts.find(p => p.slug === params.slug)

  if (!post) {
    return {
      title: '文章未找到 | 生物信息文件夹'
    }
  }

  return {
    title: `${post.title} | 生物信息文件夹`,
    description: post.excerpt,
    openGraph: {
      title: post.title,
      description: post.excerpt,
      type: 'article',
      publishedTime: post.date,
      tags: post.tags
    }
  }
}

export default async function PostPage({ params }: Props) {
  const posts = await getAllPosts()
  const post = posts.find(p => p.slug === params.slug)

  if (!post) {
    notFound()
  }

  return (
    <article className="max-w-4xl mx-auto px-4 py-12">
      <header className="mb-8">
        <h1 className="text-3xl font-bold mb-4 text-gray-900 dark:text-gray-100">
          {post.title}
        </h1>
        <div className="flex flex-wrap gap-4 text-sm text-gray-600 dark:text-gray-400">
          <time dateTime={post.date}>{post.date}</time>
          <div className="flex gap-2">
            {post.tags.map(tag => (
              <span key={tag}>#{tag}</span>
            ))}
          </div>
        </div>
      </header>
      <div 
        className="prose dark:prose-invert max-w-none"
        dangerouslySetInnerHTML={{ __html: post.content }} 
      />
    </article>
  )
} 