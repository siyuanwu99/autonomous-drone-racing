# FAST-Drone-Racing


将汪博后端引用到了 drone racing 问题上，根据门的位置修改了生成 corridor 的方式，得到了较优的结果。


## Racing tracks

这里参考宋运龙的IROS [论文](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwjcycCKkJ7yAhWQBd4KHdu5CmEQwqsBegQICBAB&url=https%3A%2F%2Fwww.youtube.com%2Fwatch%3Fv%3DHebpmadjqn8&usg=AOvVaw0TBMT96kw8zh-CeKYqTCH-)  实现了 AlphaPilot 的赛道，储存到了一个`.txt` 文件中。该文件位于`poly_planner/config/AlphaPilot.txt` 

文件格式如下所示

- 第一行存六个数据，分别为初始点和结束点的 x y z 坐标
- 第二行开始每行存4个位置，分别为门的  x y z 坐标和门的朝向。朝向使用与 x 轴的夹角度数表示。这里认为门只绕 z 轴进行旋转，不绕 x y 进行旋转。（门只转 yaw 角，不转 row 和 pitch）


```
-8     -15     0     39   12     0
5     -16    0   0
40    -14    0   50
```


## Implementations 

直接修改了汪博后端中的 `sfc_gen.cpp` 函数，按照已知门的位置生成飞行走廊。该飞行走廊分成两个部分，一部分是门附近的走廊生成，另外一部分是门之间的走廊生成。

1. 对于门附近的走廊，按照固定的长宽高生成飞行走廊。
2. 对于门之间的走廊，x y 方向上的生成办法分成如下两种情况，一种是两个门的夹角小于90度的情形，以两门为四边形的对边，两门平面各延长一定距离后再连线形成平行四边形。另外一种情况是两门的夹角大于 90 度，此时若采用上述办法会导致平行四边形过于狭窄，不利于优化出合适轨迹。此处采取的办法是以两门延长后作为平行四边形的邻边，再画出两个平行面组成平行四边形。z 方向上按照两门 z 方向上距离的一定比例生成上下两个面。

注意：**飞行走廊的每个面都由朝外的单位法向量和面上一点组成**。


## TODO
- [ ] 汪博后端的微调和改进
- [ ] 调参
- [ ] 针对 split-S 的 180度大转弯情形设计飞行走廊
- [ ] 实现其它赛道
- [ ] 在更多赛道上和宋运龙方法进行对比