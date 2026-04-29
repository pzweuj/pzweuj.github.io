'use client'

import { useState } from 'react'
import Link from 'next/link'
import { usePathname } from 'next/navigation'
import { ThemeToggle } from '../ui/ThemeToggle'
import SearchDialog from '../ui/SearchDialog'

export default function Header() {
  const pathname = usePathname()
  const [mobileMenuOpen, setMobileMenuOpen] = useState(false)

  const navigation = [
    { name: '首页', href: '/' },
    { name: '归档', href: '/archives' },
    { name: 'Schema进度', href: '/schema-progress' },
    { name: '实践', href: '/projects' },
    { name: '关于', href: '/about' }
  ]

  return (
    <header className="sticky top-0 z-50 w-full border-b border-gray-200 dark:border-gray-800 bg-white/80 dark:bg-gray-900/80 backdrop-blur-sm">
      <nav className="max-w-4xl mx-auto px-4 flex justify-between items-center h-14 md:h-16">
        {/* 左侧：移动端汉堡 + 桌面端导航 */}
        <div className="flex items-center">
          {/* 移动端汉堡按钮 - 最小 44px 触摸区域 */}
          <button
            className="md:hidden flex items-center justify-center w-11 h-11 -ml-2 text-gray-600 dark:text-gray-300 hover:text-gray-900 dark:hover:text-white rounded-lg transition-colors"
            onClick={() => setMobileMenuOpen(!mobileMenuOpen)}
            aria-label={mobileMenuOpen ? '关闭菜单' : '打开菜单'}
            aria-expanded={mobileMenuOpen}
          >
            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              {mobileMenuOpen ? (
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
              ) : (
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 6h16M4 12h16M4 18h16" />
              )}
            </svg>
          </button>

          {/* 桌面端导航 */}
          <div className="hidden md:flex items-center space-x-6">
            {navigation.map((item) => (
              <Link
                key={item.name}
                href={item.href}
                className={`text-sm font-medium transition-colors hover:text-blue-600 dark:hover:text-blue-400 ${
                  pathname === item.href
                    ? 'text-blue-600 dark:text-blue-400'
                    : 'text-gray-600 dark:text-gray-300'
                }`}
              >
                {item.name}
              </Link>
            ))}
          </div>
        </div>

        {/* 右侧：搜索 + 主题切换 */}
        <div className="flex items-center gap-1">
          <SearchDialog />
          <ThemeToggle />
        </div>
      </nav>

      {/* 移动端菜单 - 带过渡动画 */}
      <div
        className={`md:hidden overflow-hidden transition-all duration-200 ease-in-out ${
          mobileMenuOpen ? 'max-h-80 opacity-100' : 'max-h-0 opacity-0'
        }`}
      >
        <div className="border-t border-gray-200 dark:border-gray-800 bg-white dark:bg-gray-900 px-4 py-2">
          {navigation.map((item) => (
            <Link
              key={item.name}
              href={item.href}
              onClick={() => setMobileMenuOpen(false)}
              className={`flex items-center px-3 py-3 rounded-lg text-base font-medium transition-colors ${
                pathname === item.href
                  ? 'text-blue-600 dark:text-blue-400 bg-blue-50 dark:bg-blue-900/20'
                  : 'text-gray-600 dark:text-gray-300 hover:text-blue-600 dark:hover:text-blue-400 hover:bg-gray-50 dark:hover:bg-gray-800'
              }`}
            >
              {item.name}
            </Link>
          ))}
        </div>
      </div>
    </header>
  )
}
