'use client'

import Link from 'next/link'
import { usePathname } from 'next/navigation'
import { ThemeToggle } from '../ui/ThemeToggle'
import SearchDialog from '../ui/SearchDialog'

export default function Header() {
  const pathname = usePathname()

  const navigation = [
    { name: '首页', href: '/' },
    { name: '归档', href: '/archives' },
    { name: '实践', href: '/projects' },
    { name: '关于', href: '/about' }
  ]

  return (
    <header className="sticky top-0 z-50 w-full border-b border-gray-200 dark:border-gray-800 bg-white dark:bg-gray-900">
      <nav className="max-w-4xl mx-auto px-4 flex justify-between items-center h-16">
        <div className="flex items-center space-x-8">
          {navigation.map((item) => (
            <Link
              key={item.name}
              href={item.href}
              className={`text-sm font-medium transition-colors ${
                pathname === item.href
                  ? 'text-blue-600 dark:text-blue-400'
                  : 'text-gray-600 dark:text-gray-300 hover:text-blue-600 dark:hover:text-blue-400'
              }`}
            >
              {item.name}
            </Link>
          ))}
        </div>
        <div className="flex items-center space-x-4">
          <SearchDialog />
          <ThemeToggle />
        </div>
      </nav>
    </header>
  )
}