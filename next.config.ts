/** @type {import('next').NextConfig} */
const nextConfig = {
  output: 'export',  // 添加这行
  basePath: '/pzweuj.github.io',
  images: {
    unoptimized: true,
  },
  eslint: {
    ignoreDuringBuilds: true,
  },
  typescript: {
    ignoreBuildErrors: true,
  },
}

module.exports = nextConfig