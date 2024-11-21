import type { NextConfig } from 'next'

const nextConfig: NextConfig = {
  output: "export",  // 启用静态导出
  reactStrictMode: true,
  images: {
    unoptimized: true, // 静态导出需要
  }
}

export default nextConfig 