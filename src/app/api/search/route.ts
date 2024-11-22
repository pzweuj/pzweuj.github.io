import { generateSearchIndex } from '@/lib/search'
import { NextResponse } from 'next/server'

export const dynamic = 'force-static'

export async function GET() {
  const searchIndex = await generateSearchIndex()
  return NextResponse.json(searchIndex)
}