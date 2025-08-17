---
layout: post
read_time: true
show_date: true
title:  使用云服务器+frp访问内网机器
date:   2025-08-17 13:40:20 -0600
header-img: img/20250114/gondolin.jpg
tags: [计算生物学,服务器连接]
author: 孙睿
mathjax: yes
catalog: true
--- 

这篇blog介绍怎样使用公网云服务器，通过frp服务实现对内网服务器的远程访问。

最近在外地研究所工作，需要访问在北京研究所的一台内网机器。尝试了远程桌面、zerotier内网穿透等多种方式后，最后决定使用腾讯云 + frp构建tunneling来实现远程访问。这个策略几乎不受两个研究所的远程桌面防火墙的限制，延迟显著低于zerotier，费用上也比自费搞一个高速的无线wifi再远程桌面连接便宜许多。只要两个研究所内网能够访问公网的腾讯云服务器即可。

基本框架：

```
笔记本（外网） → 腾讯云服务器（frp server） ← frp client（研究所服务器）
                                               ↖ SSH
```

### 云服务器配置
首先注册一个腾讯云服务器。登陆后做如下设置

```bash
# 示例（以 Linux amd64 为例，查看最新版本：https://github.com/fatedier/frp/releases）
wget https://github.com/fatedier/frp/releases/latest/download/frp_0.52.3_linux_amd64.tar.gz
tar -zxpf frp_0.52.3_linux_amd64.tar.gz
cd frp_0.52.3_linux_amd64
``` 

之后配置frps.ini文件

```bash
[common]
bind_port = 443  # FRP 服务端端口, 
#这个端口要考虑两个研究所的网络防火墙限制，443这种常用端口通常ok
token = your_secure_token  # 访问认证 token，非常重要，请自定义一个复杂字符串
``` 

之后运行frps服务

```bash
nohup ./frps -c frps.ini &
```

### 研究所服务器配置 

同样在研究所机器上下载frp，之后配置frpc.ini文件

```bash
[common]
server_addr = 腾讯云服务器公网IP
server_port = 443 # 和服务端的 bind_port 一致
token = your_secure_token # 和服务端的 token 一致

[ssh]
type = tcp
local_ip = ********* # 研究所服务器的local ip
local_port = 22 # 内网服务器的SSH端口，通常为22
remote_port = 6000 # 映射到腾讯云服务器的端口，自定义
```

之后运行frpc 

```bash
nohup ./frpc -c frpc.ini &
```

### ssh连接 

现在可以使用ssh连接到北京研究所内的服务器

```bash
ssh -p 6000 username@腾讯云服务器公网IP

# username 为北京服务器的用户名
```

为了方便在vscode中直接ssh连接，修改本地ssh配置。打开VSCODE, ctrl + shift + p 打开remote-ssh configure，在其中添加

```bash
# 北京实验室服务器 (通过腾讯云FRP隧道)
Host A800-tencent
    HostName 腾讯云公网ip
    User 北京实验室用户名
    Port 6000
    IdentityFile ~/.ssh/id_rsa
```

修改之后可以直接在vscode里面remote-ssh连接。

如果需要配置ssh密钥免密码登录。
```bash
mkdir -p ~/.ssh
chmod 700 ~/.ssh

vim ~/.ssh/authorized_keys
chmod 600 ~/.ssh/authorized_keys
```
在authorized_keys中添加笔记本电脑的公钥~/.ssh/id_rsa.pub中的公钥，每个公钥占一行。

以上就是使用腾讯云远程连接内网机器的策略。感谢北大的逸鸣师兄提供frp解决方案思路，以及gemini2.5-flash提供的实现流程。