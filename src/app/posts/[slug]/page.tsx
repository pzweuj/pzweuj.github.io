import { getAllPosts } from '@/lib/markdown'
import { notFound } from 'next/navigation'
import type { Metadata } from 'next'
import selfConfig from '@/config/self.config'
import 'katex/dist/katex.min.css'

interface Props {
  params: Promise<{ slug: string }>
}

export async function generateMetadata({ params }: Props): Promise<Metadata> {
  const { slug } = await params
  const posts = await getAllPosts()
  const post = posts.find(p => p.slug === slug)

  if (!post) {
    return {
      title: '文章未找到'
    }
  }

  return {
    title: post.title,
    description: post.excerpt,
    openGraph: {
      title: post.title,
      description: post.excerpt,
      type: 'article',
      publishedTime: post.date,
      tags: post.tags,
      url: `${selfConfig.siteUrl}/posts/${post.slug}`,
    },
    twitter: {
      card: 'summary',
      title: post.title,
      description: post.excerpt,
    },
    alternates: {
      canonical: `${selfConfig.siteUrl}/posts/${post.slug}`,
    },
  }
}

export default async function PostPage({ params }: Props) {
  const { slug } = await params
  const posts = await getAllPosts()
  const post = posts.find(p => p.slug === slug)

  if (!post) {
    notFound()
  }

  const jsonLd = {
    '@context': 'https://schema.org',
    '@type': 'BlogPosting',
    headline: post.title,
    datePublished: post.date,
    description: post.excerpt,
    keywords: post.tags,
    url: `${selfConfig.siteUrl}/posts/${post.slug}`,
    author: {
      '@type': 'Person',
      name: selfConfig.author,
      url: selfConfig.social.github,
    },
    publisher: {
      '@type': 'Person',
      name: selfConfig.author,
    },
    mainEntityOfPage: {
      '@type': 'WebPage',
      '@id': selfConfig.siteUrl,
    },
  }

  return (
    <article className="max-w-4xl mx-auto px-4 py-12">
      <script
        type="application/ld+json"
        dangerouslySetInnerHTML={{ __html: JSON.stringify(jsonLd) }}
      />
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
        className="prose max-w-none"
        dangerouslySetInnerHTML={{ __html: post.content }}
      />
    </article>
  )
}

// 添加 generateStaticParams 函数来生成所有可能的文章路径
export async function generateStaticParams() {
  const posts = await getAllPosts()
  return posts.map((post) => ({
    slug: post.slug
  }))
} 