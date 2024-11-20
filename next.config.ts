/** @type {import('next').NextConfig} */
const nextConfig = {
  output: 'export',  // 添加这行
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
