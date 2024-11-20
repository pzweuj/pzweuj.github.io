import { getAllPosts } from '@/lib/markdown'
import { NextResponse } from 'next/server'

export const dynamic = 'force-static'
export const revalidate = false

// 预生成包含所有文章的响应
export async function GET() {
  const posts = await getAllPosts()
  const searchData = posts.map(post => ({
    title: post.title,
    excerpt: post.excerpt,
    slug: post.slug,
    date: post.date
  }))
  
  return NextResponse.json(searchData)
} 