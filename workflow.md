# Workflow
# 摘要
本文件简单叙述了工作流，详细的操作可参考[GitHub pull request入门（图解+原理+git命令+可能有用的经验） - 知乎](https://zhuanlan.zhihu.com/p/672447698)
## 零、Git 安装和配置以及 Github 创建账号
可上网查询，此处省略
## 一、 初始化
访问 [WinstonMeursault/SDU-ASC26: ASC26 repo of SDU](https://github.com/WinstonMeursault/SDU-ASC26) 并点击右上角的 Fork 按钮。之后访问自己仓库的 SDU-ASC26 仓库，点击右上角绿色的 Code 按钮，在下拉框中复制 HTTPS 标签页下的那一串网址，之后在本地某个文件夹（作为项目文件夹的上一级）处运行命令：
```bash
git clone 网址
```
这个命令会新建一个`SDU-ASC26`目录并将仓库中的所有文件克隆到这个本地目录中。
之后执行
```bash
git remote add upstream https://github.com/WinstonMeursault/SDU-ASC26
```
以建立与原项目的连接。

# 二、对项目进行修改
执行
```bash
git switch main
git fetch upstream
git merge upstream/main
git push origin main
git switch -c 你的分支名
```
新建分支（分支名可自己起，如`dev`）。
~~然后就没什么好说的，修改项目~~
~~推荐用[vscode](https://code.visualstudio.com/)，因为我也在用（）~~

# 三、提交修改
运行
```bash
git add .
git commit -m "commit信息，最好把自己的修改大致地写进来"
git fetch upstream main
git rebase upstream/main
```
如果有代码冲突，先人工解决掉，再执行上面第二个命令。
接下来运行
```bash
git push origin 你的分支名
```
提交更改。
# 四、提 Pull Request
进入 Github 里面自己 Fork 的项目，会发现绿色按钮下面多出了一行字，点击 Contribute 中的 Open Pull Request ，将自己的分支提交到远程的`main`分支即可。
# 五、删除自己的分支
在 Pull Request 接受后，运行
```bash
git switch main
git fetch upstream
git merge upstream/main
git push origin main
git branch -d 你的分支名（不是 main）
git push origin --delete 你的分支名（不是 main）
```
删除自己的分支，之后重复第二至五步。
# 参考
1. [GitHub pull request入门（图解+原理+git命令+可能有用的经验） - 知乎](https://zhuanlan.zhihu.com/p/672447698)
2. [十分钟学会正确的github工作流，和开源作者们使用同一套流程](https://www.bilibili.com/video/BV19e4y1q7JJ?vd_source=90d8ee5c476dd6c93bfe95eba=aea8302)