# LumpedMassMethodMooring
Use the lump mass method to calculate cable forces

使用集中质量法求解悬链线的缆绳力

算法很简单，主要思路：
1. 由锚固点坐标推导出锚链刚好被完全拉起时的拉力与系泊点坐标的关系
2. 根据悬链线长度 l<sub>p</sub>、两点直线长度 l<sub>s</sub> 与锚链长度 l<sub>m</sub> 的关系可分为三种情况求解力
    - l<sub>m</sub>>l<sub>p</sub>:此时锚链处于拖地状态
    - l<sub>p</sub>>l<sub>m</sub>>l<sub>s</sub>:锚链完全拉起，无拖地长度
    - l<sub>s</sub>>l<sub>m</sub>:锚链几乎伸直，可直接按胡克定律求解

## 更新记录

### 20200301

#### 变更
* 从论文中提取源代码并在OpenFOAM-3.0.1中编译。
* **代码来自`于含`的硕士学位论文**[系泊载液浮体水动力特性的数值及试验研究](http://gb.oversea.cnki.net/KCMS/detail/detailall.aspx?filename=1018869021.nh&dbcode=CMFD&dbname=CMFD2019)
* 基于`GNU`许可，开源

### 20200302

#### 内容

* 对论文中验证算例的验证



### 20200303-20200420

#### 咕咕咕

* 忙其他东西去了

### 20200421

#### 优化

* 将代码各行进行解析注释 
* 同时开始思考算法如何扩展到三维

### 20200423

#### 新增

* 将算法扩展到了三维 
* 发布了编译成功的第一版

### 20200427

#### 修复

* 对之前的二维验证算例再度进行计算，发现巨大bug
* 现已修复，并发布

