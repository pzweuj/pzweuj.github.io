import { getProjectDocs, getProjectIndex } from '@/lib/projects'
import { ProjectSidebar } from '@/components/projects/ProjectSidebar'

export default async function ProjectsPage() {
  const chapters = await getProjectDocs()
  const index = await getProjectIndex()
  
  return (
    <div className="flex">
      {/* 侧边栏 */}
      <ProjectSidebar chapters={chapters} />
      
      {/* 主内容区 */}
      <main className="flex-1 max-w-4xl mx-auto px-4 py-12">
        <div className="prose dark:prose-invert max-w-none">
          <div dangerouslySetInnerHTML={{ __html: index.content }} />
        </div>
      </main>
    </div>
  )
} 