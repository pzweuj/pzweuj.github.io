---
title: Tauri2文件选择器爬坑
tags: software
---

我想实现的功能是把项目封装后，点击按钮保存文件时打开一个文件对话框，然后可以自行选择保存路径，设定文件名称，然后保存。

使用的是下面这个代码

```tsx
import { save } from '@tauri-apps/plugin-dialog'

const filePath = await save({
    defaultPath: `maneloca_${version}_${new Date().toISOString().slice(0, 10)}.tsv`, // 默认文件名
    filters: [{
      name: 'TSV Files',
      extensions: ['tsv'] // 限制文件类型为 TSV
    }]
})
```

但是封装后点击对应的按钮，没有反应。因为我没有任何前端知识，和AI斗智斗勇了N个回合，还是无法找到合适的方案。后来发现应该是AI的知识还停留在Tauri 1，只能自己查文档了。

## 保存至文件对话框

下面是Tauri 2保存文件的[示例代码](https://tauri.org.cn/plugin/dialog/#save-to-file-dialog)。

```tsx
import { save } from '@tauri-apps/plugin-dialog';
// when using `"withGlobalTauri": true`, you may use
// const { save } = window.__TAURI__.dialog;

// Prompt to save a 'My Filter' with extension .png or .jpeg
const path = await save({
  filters: [
    {
      name: 'My Filter',
      extensions: ['png', 'jpeg'],
    },
  ],
});
console.log(path);
// Prints the chosen path
```

但是呢，我的代码本就是按照这个方案进行的，仍然无法弹出对应的对话框。在对代码补充了一个try之后，发现这段代码应该是报错了，没有跑进去，然而后端没有错误提示。AI说是tauri配置文件的问题。

## 配置文件

AI说要先在配置文件中允许插件运行，给我修改建议如下所示。

```json
{
  "tauri": {
    "allowlist": {
      "dialog": {
        "all": true, // 允许所有 dialog 功能
        "save": true // 明确允许 save 功能
      }
    }
  }
}
```

但我在修改时，提示配置文件不允许"tauri"字段，感觉又是版本迭代后配置的方式不同了。在[Tauri 2的文档](https://v2.tauri.app/zh-cn/develop/plugins/)中写了在需要使用插件的项目中应该如何配置。

我用了tauri的fs和dialog插件，那么就应该在配置文件中加入

```json
"plugins": {
    "dialog" : true,
    "fs": true
}
```

我这样设置后，还是不行。

## 需要安装

之前我一直是只安装了js库

```cmd
npm install @tauri-apps/plugin-fs
```

实际上要通过tauri安装

```cmd
npm run tauri add dialog
```

所以其实就这样然后其他配置都不用改，AI误我😭。


最后测试可以用了

```cmd
npm run tauri dev
```

我觉得有点复杂，考虑转投electron。
