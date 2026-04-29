import { getAllSchemaProgress } from '@/lib/markdown'
import { notFound } from 'next/navigation'
import { Metadata } from 'next'
import 'katex/dist/katex.min.css'

interface Props {
  params: Promise<{ slug: string }>
}

export async function generateMetadata({ params }: Props): Promise<Metadata> {
  const { slug } = await params
  const posts = await getAllSchemaProgress()
  const post = posts.find(p => p.slug === slug)

  if (!post) {
    return {
      title: '进度未找到 | Schema进度'
    }
  }

  return {
    title: `${post.title} | Schema进度`,
    description: post.excerpt,
  }
}

export default async function SchemaProgressPostPage({ params }: Props) {
  const { slug } = await params
  const posts = await getAllSchemaProgress()
  const post = posts.find(p => p.slug === slug)

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
        </div>
      </header>
      <div
        className="prose dark:prose-invert max-w-none"
        dangerouslySetInnerHTML={{ __html: post.content }}
      />
    </article>
  )
}

export async function generateStaticParams() {
  const posts = await getAllSchemaProgress()
  return posts.map((post) => ({
    slug: post.slug
  }))
}