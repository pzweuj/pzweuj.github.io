import { GithubIcon } from '../ui/Icons'
import selfConfig from '@/config/self.config'

// 添加 EmailIcon 组件
const EmailIcon = () => (
  <svg
    xmlns="http://www.w3.org/2000/svg"
    width="24"
    height="24"
    viewBox="0 0 24 24"
    fill="none"
    stroke="currentColor"
    strokeWidth="2"
    strokeLinecap="round"
    strokeLinejoin="round"
  >
    <rect x="2" y="4" width="20" height="16" rx="2"></rect>
    <path d="m22 7-8.97 5.7a1.94 1.94 0 0 1-2.06 0L2 7"></path>
  </svg>
)

export default function Footer() {
  return (
    <footer className="border-t bg-gray-50 dark:bg-gray-900 dark:border-gray-800">
      <div className="max-w-6xl mx-auto px-4 py-8">
        <div className="flex flex-col items-center justify-center space-y-4">
          {/* 社交链接 */}
          <div className="flex space-x-6">
            <a
              href={selfConfig.social.github}
              target="_blank"
              rel="noopener noreferrer"
              className="text-gray-600 hover:text-gray-900 dark:text-gray-400 dark:hover:text-white"
              aria-label="GitHub"
            >
              <GithubIcon />
            </a>
            <a
              href={`mailto:${selfConfig.social.email}`}
              className="text-gray-600 hover:text-gray-900 dark:text-gray-400 dark:hover:text-white"
              aria-label="Email"
            >
              <EmailIcon />
            </a>
          </div>

          {/* 版权信息 */}
          <div className="text-center text-sm text-gray-600 dark:text-gray-400">
            <p>© {new Date().getFullYear()} {selfConfig.author}. All rights reserved.</p>
            <p className="mt-1">
              Built with{' '}
              <a 
                href="https://nextjs.org"
                target="_blank"
                rel="noopener noreferrer"
                className="text-blue-600 hover:text-blue-500 dark:text-blue-400"
              >
                Next.js
              </a>
            </p>
          </div>
        </div>
      </div>
    </footer>
  )
} 