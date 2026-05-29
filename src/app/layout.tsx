import { Inter } from 'next/font/google'
import type { Metadata } from 'next'
import { ThemeProvider } from '@/components/providers/ThemeProvider'
import Header from '@/components/layout/Header'
import Footer from '@/components/layout/Footer'
import ScrollToTop from '@/components/ui/ScrollToTop'
import { CodeCopyScript } from '@/components/ui/CodeCopyScript'
import selfConfig from '@/config/self.config'
import './globals.css'

const inter = Inter({ subsets: ['latin'] })

export const metadata: Metadata = {
  title: {
    default: selfConfig.title,
    template: `%s | ${selfConfig.title}`,
  },
  description: selfConfig.description,
  authors: [{ name: selfConfig.author, url: selfConfig.social.github }],
  creator: selfConfig.author,
  metadataBase: new URL(selfConfig.siteUrl),
  icons: {
    icon: '/icon.ico',
    shortcut: '/icon.ico',
  },
  referrer: 'strict-origin-when-cross-origin',
  formatDetection: {
    email: false,
    address: false,
    telephone: false,
  },
  openGraph: {
    type: 'website',
    locale: 'zh_CN',
    url: selfConfig.siteUrl,
    siteName: selfConfig.title,
    title: selfConfig.title,
    description: selfConfig.description,
  },
  twitter: {
    card: 'summary',
    title: selfConfig.title,
    description: selfConfig.description,
  },
  alternates: {
    canonical: selfConfig.siteUrl,
    types: {
      'application/rss+xml': [
        { url: '/rss.xml', title: 'RSS Feed' },
      ],
      'application/atom+xml': [
        { url: '/atom.xml', title: 'Atom Feed' },
      ],
    },
  },
}

export default function RootLayout({
  children,
}: {
  children: React.ReactNode
}) {
  return (
    <html lang={selfConfig.language} suppressHydrationWarning>
      <head>
        <meta name="theme-color" content="#ffffff" media="(prefers-color-scheme: light)" />
        <meta name="theme-color" content="#111827" media="(prefers-color-scheme: dark)" />
        <script
          type="application/ld+json"
          dangerouslySetInnerHTML={{
            __html: JSON.stringify({
              '@context': 'https://schema.org',
              '@type': 'WebSite',
              name: selfConfig.title,
              description: selfConfig.description,
              url: selfConfig.siteUrl,
              author: {
                '@type': 'Person',
                name: selfConfig.author,
                url: selfConfig.social.github,
              },
              inLanguage: selfConfig.language,
            }),
          }}
        />
      </head>
      <body className={`${inter.className} min-h-screen flex flex-col bg-white dark:bg-gray-900`}>
        <ThemeProvider
          attribute="class"
          defaultTheme="system"
          enableSystem
          disableTransitionOnChange
        >
          <Header />
          <main className="flex-1">
            {children}
          </main>
          <Footer />
          <ScrollToTop />
          <CodeCopyScript />
        </ThemeProvider>
      </body>
    </html>
  )
}
