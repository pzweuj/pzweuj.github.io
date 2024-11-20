'use client'

import { useState, useEffect } from 'react'
import Link from 'next/link'
import { usePathname } from 'next/navigation'
import type { ProjectChapter } from '@/lib/projects'
import { Menu, X } from 'lucide-react'

interface ProjectSidebarProps {
  chapters: ProjectChapter[]
}

export function ProjectSidebar({ chapters }: ProjectSidebarProps) {
  const pathname = usePathname()
  const [expandedChapters, setExpandedChapters] = useState<string[]>([])
  const [isMobileMenuOpen, setIsMobileMenuOpen] = useState(false)
  
  // 从 localStorage 加载展开状态
  useEffect(() => {
    const saved = localStorage.getItem('expandedChapters')
    if (saved) {
      setExpandedChapters(JSON.parse(saved))
    }
  }, [])
  
  // 保存展开状态到 localStorage
  const toggleChapter = (chapterId: string) => {
    const newExpanded = expandedChapters.includes(chapterId)
      ? expandedChapters.filter(id => id !== chapterId)
      : [...expandedChapters, chapterId]
    
    setExpandedChapters(newExpanded)
    localStorage.setItem('expandedChapters', JSON.stringify(newExpanded))
  }
  
  // 导航栏内容组件
  const SidebarContent = () => (
    <nav className="p-4 space-y-8">
      <div>
        <h3 className="font-medium text-gray-900 dark:text-gray-100 mb-2">
          <Link 
            href="/projects"
            className={`inline-flex items-center space-x-3 ${
              pathname === '/projects'
                ? 'text-blue-600 dark:text-blue-400'
                : 'text-gray-600 dark:text-gray-300 hover:text-blue-600 dark:hover:text-blue-400'
            }`}
          >
            <svg 
              className="w-5 h-5" 
              fill="none" 
              viewBox="0 0 24 24" 
              stroke="currentColor"
            >
              <path 
                strokeLinecap="round" 
                strokeLinejoin="round" 
                strokeWidth={2} 
                d="M3 12l2-2m0 0l7-7 7 7M5 10v10a1 1 0 001 1h3m10-11l2 2m-2-2v10a1 1 0 01-1 1h-3m-6 0a1 1 0 001-1v-4a1 1 0 011-1h2a1 1 0 011 1v4a1 1 0 001 1m-6 0h6" 
              />
            </svg>
            <span className="text-base">返回首页</span>
          </Link>
        </h3>
      </div>
      
      {chapters.map(chapter => (
        <div key={chapter.id}>
          <button
            onClick={() => toggleChapter(chapter.id)}
            className="w-full flex items-center justify-between font-medium text-gray-900 dark:text-gray-100 mb-2 group hover:text-blue-600 dark:hover:text-blue-400"
          >
            <span>{chapter.title}</span>
            <svg
              className={`w-4 h-4 transition-transform duration-200 ${
                expandedChapters.includes(chapter.id) ? 'transform rotate-180' : ''
              }`}
              fill="none"
              viewBox="0 0 24 24"
              stroke="currentColor"
            >
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                strokeWidth={2}
                d="M19 9l-7 7-7-7"
              />
            </svg>
          </button>
          {expandedChapters.includes(chapter.id) && (
            <ul className="space-y-2 ml-4">
              {chapter.docs.map(doc => (
                <li key={doc.slug}>
                  <Link
                    href={`/projects/${doc.slug}`}
                    className={`block text-sm ${
                      pathname === `/projects/${doc.slug}`
                        ? 'text-blue-600 dark:text-blue-400'
                        : 'text-gray-600 dark:text-gray-300 hover:text-blue-600 dark:hover:text-blue-400'
                    }`}
                    onClick={() => setIsMobileMenuOpen(false)}
                  >
                    {doc.title}
                  </Link>
                </li>
              ))}
            </ul>
          )}
        </div>
      ))}
    </nav>
  )
  
  return (
    <>
      {/* 桌面端侧边栏 */}
      <aside className="hidden md:block w-64 h-[calc(100vh-4rem)] overflow-y-auto sticky top-16 border-r border-gray-200 dark:border-gray-800">
        <SidebarContent />
      </aside>
      
      {/* 移动端汉堡菜单按钮 */}
      <button
        onClick={() => setIsMobileMenuOpen(!isMobileMenuOpen)}
        className="md:hidden fixed bottom-6 left-6 z-50 p-3 bg-blue-600 text-white rounded-full shadow-lg hover:bg-blue-700 transition-colors"
      >
        {isMobileMenuOpen ? <X size={24} /> : <Menu size={24} />}
      </button>
      
      {/* 移动端侧边栏 */}
      {isMobileMenuOpen && (
        <div className="md:hidden fixed inset-0 z-40 bg-white dark:bg-gray-900">
          <div className="h-full overflow-y-auto">
            <SidebarContent />
          </div>
        </div>
      )}
    </>
  )
} 