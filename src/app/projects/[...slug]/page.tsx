import { getProjectDocs } from '@/lib/projects'
import { ProjectSidebar } from '@/components/projects/ProjectSidebar'
import { notFound } from 'next/navigation'

interface Props {
  params: {
    slug: string[]
  }
}

// 添加这个函数来生成静态路径
export async function generateStaticParams() {
  const chapters = await getProjectDocs()
  
  // 获取所有文档的 slug
  const paths = chapters.flatMap(chapter => 
    chapter.docs.map(doc => ({
      slug: doc.slug.split('/')  // 将 slug 字符串分割成数组
    }))
  )
  
  return paths
}

export default async function ProjectDocPage({ params }: Props) {
  const chapters = await getProjectDocs()
  const slug = params.slug.join('/')
  
  // 查找当前文档
  const doc = chapters
    .flatMap(chapter => chapter.docs)
    .find(doc => doc.slug === slug)
  
  if (!doc) {
    notFound()
  }
  
  return (
    <div className="flex">
      <ProjectSidebar chapters={chapters} />
      
      <main className="flex-1 max-w-4xl mx-auto px-4 py-12">
        <div className="prose dark:prose-invert max-w-none">
          <h1>{doc.title}</h1>
          <div dangerouslySetInnerHTML={{ __html: doc.content }} />
        </div>
      </main>
    </div>
  )
} 