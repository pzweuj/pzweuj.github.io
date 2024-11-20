import { getAllPosts } from '@/lib/markdown'
import { NextResponse } from 'next/server'

export async function GET(request: Request) {
  const { searchParams } = new URL(request.url)
  const query = searchParams.get('q')?.toLowerCase()

  if (!query) {
    return NextResponse.json([])
  }

  const posts = await getAllPosts()
  
  const results = posts
    .filter(post => post.title.toLowerCase().includes(query))
    .map(post => ({
      title: post.title,
      excerpt: post.excerpt,
      slug: post.slug,
      date: post.date
    }))

  return NextResponse.json(results)
} 