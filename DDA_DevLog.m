% 此文件为开发日志，记录各种改动
% 20170405 完成了DDA_Compute 、 DDA_TargetOrienGener
% 
% DDA_TargetOrienGener 可以均匀输出球面上的取向点，以及可以附带颗粒物的内在旋转
% DDA_Compute 按照对应的颗粒物取向，输出颗粒物的P阵（P为极化向量）

%0409 DDA_ResultOutput的检查、优化
%目前无问题 以下均为对 2586个探测点的计算结果
%最初版本，需要25860s左右的时间计算出结果
%优化部分
%1 做了vertiE 和 horiE一起在循环内部计算的优化
%2 生成各种temp变量 ， 极大的减少了deg2rad squeeze 等函数的调用次数 
%此时耗时接近3000s
%进行parfor优化
%做了parfor优化 现在输出一次结果耗时1000s左右 （cpu i5 6500） 4核机器下，速度提升至运来的3倍左右，因为存在通信开销

%0411 DDA_SphereFigure 函数完成 用于显示出射muller矩阵元的球面谱
%DDA_SphereFigure 函数优化，保证了相同大小的元素显示的颜色为同一的

%0412 DDA_ResultOutput的优化
%采用矩阵计算、矩阵乘法的办法进行优化
%编写了两个新的结果输出脚本 DDA_ResultOutputMatrixVersion 和 DDA_ResultOutputParforMatrixVersion
%此两脚本均采用了矩阵运算的方法来加速计算
%两种计算结果所需时间
%DDA_ResultOutputMatrixVersion 14.0s
%DDA_ResultOutputParforMatrixVersion 12.5s
%经过观察发现，尽管没有parfor，DDA_ResultOutputMatrixVersion会自动调用多核，使得其效率与parfor版本一样
%其输出结果与DDA_ResultOutput有10e-6的误差，造成此原因的结果可能是DDA_ResultOutput在for循环中的round off现象

%0413 ResultOutput 加入了两种GPU算法 ResultOutputGPUVersion 和 ResultOutputParforGPUVersion 
%在面对大数据时，ResultOutputGPUVersion 运行时间为 DDA_ResultOutputMatrixVersion 的三分之一
%ResultOutputParforGPUVersion 运行时间为 DDA_ResultOutputMatrixVersion 的十二分之一
%以上脚本均以验证， 其输出结果与DDA_ResultOutput一样

% 以DDA_ResultOutput的运行效率为1 各结果输出脚本运行效率如下
% DDA_ResultOutputMatrixVersion = DDA_ResultOutputParforMatrixVersion  = 80
% ResultOutputGPUVersion = 240
% ResultOutputParforGPUVersion = 960

%0416 使用了球的mie散射来验证整个DDA的正确性
%球的散射矩阵，在球面谱上 S11 S22 S23 S32 S33 S44为主要值
%主要值上，mie散射的结果和DDA的结果是类似的 
%deltaPercent = MatrixCompare( mieTempMuller(:,indice), tempMuller(:,indice) ) 命令
%返回的结果显示，deltaPercent一般相差为0.2左右，这说明程序的结果是对的，自己撰写的DDA是没有问题的。
%但是0.2也是不小的数字，这说明这个DDA的结果可能与真实情况有一定的相差。

%0417 开发完了DDA_MullerT_PhiFillUp 已用球散射验证其正确性
%在2586个探测点的情况下，运行旋转程序导致的偏差为 10e-8
%另外，对于粒子的旋转t_phi，对于探测角的旋转L_phi，都是以激光入射方向（作为轴矢量）为正方向

%修改了DDA_Compute中 旋转部分deg与rad的误用情况 跑了球的有取向有内转旋转的结果 初步认为整个程序没有问题

%0419 开始修改DDA_Compute DDA_ConvAccelerate  核心是把 Ag 直接先存储在显存里 避免将A从内存传入至显存中
%这种修改之后，DDA_ConvAccelerate  的接口从接受Af 变为接受Afg 故而DDA_Compute需要做对应的修改
%将outArrayg 和 Yg 写死在显存里，此函数效率得到了2.5倍的提升
%修改后的函数和原来有e-4的差别，但是用球散射验证了其正确性

%将 A的生成放在 eDirection 之外 ， 节省生成A的时间 ，且逻辑更清晰

% 0420 开发了DDA_MemoryEstimate，能用来估算需要消耗的内存、显存
% 测试了 4.5 2.2 2.2um粒子的计算速度 ， 本模拟程序是DDSCAT的10倍以上
% 重写了 DDA_ModelGenerate ，现在此函数的功能是直接输入粒子的几何线度，然后自适应生成点阵。

% 0421
% 重写了 DDA_StrucPlot 此函数会在两个版本的ModelGenerate函数中调用，起到生成玩Model之后，将粒子的
% 空间几何图像绘出的作用

% 0424
% 完成了 DDA_PointMuller函数，可以通过拟合生成球面上任意点的muller矩阵