import type { MetadataRoute } from 'next'
import { getAllPosts, getAllSchemaProgress } from '@/lib/markdown'
import { getProjectDocs } from '@/lib/projects'
import selfConfig from '@/config/self.config'

export const dynamic = 'force-static'

export default async function sitemap(): Promise<MetadataRoute.Sitemap> {
  const posts = await getAllPosts()
  const chapters = await getProjectDocs()
  const schemaProgress = await getAllSchemaProgress()

  const postEntries: MetadataRoute.Sitemap = posts.map((post) => ({
    url: `${selfConfig.siteUrl}/posts/${post.slug}`,
    lastModified: new Date(post.date),
    changeFrequency: 'monthly' as const,
    priority: 0.8,
  }))

  const projectEntries: MetadataRoute.Sitemap = chapters.flatMap((chapter) =>
    chapter.docs.map((doc) => ({
      url: `${selfConfig.siteUrl}/projects/${doc.slug}`,
      lastModified: new Date(),
      changeFrequency: 'monthly' as const,
      priority: 0.6,
    }))
  )

  const schemaEntries: MetadataRoute.Sitemap = schemaProgress.map((post) => ({
    url: `${selfConfig.siteUrl}/schema-progress/${post.slug}`,
    lastModified: new Date(post.date),
    changeFrequency: 'weekly' as const,
    priority: 0.6,
  }))

  return [
    {
      url: selfConfig.siteUrl,
      lastModified: new Date(),
      changeFrequency: 'daily',
      priority: 1,
    },
    {
      url: `${selfConfig.siteUrl}/archives`,
      lastModified: new Date(),
      changeFrequency: 'daily',
      priority: 0.7,
    },
    {
      url: `${selfConfig.siteUrl}/projects`,
      lastModified: new Date(),
      changeFrequency: 'monthly',
      priority: 0.7,
    },
    {
      url: `${selfConfig.siteUrl}/schema-progress`,
      lastModified: new Date(),
      changeFrequency: 'weekly',
      priority: 0.7,
    },
    {
      url: `${selfConfig.siteUrl}/about`,
      lastModified: new Date(),
      changeFrequency: 'monthly',
      priority: 0.5,
    },
    ...postEntries,
    ...projectEntries,
    ...schemaEntries,
  ]
}
