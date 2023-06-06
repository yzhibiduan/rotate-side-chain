# rotate-side-chain
这是一个对单个残基旋转侧链二面角的简单脚本，旋转后可能会发生原子碰撞等问题，请仔细使用。
运行环境python==3.9 biopython==1.81 numpy==1.21.2，其他环境未经过测试。
二面角定义来自于http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
注意！此脚本可能无法正确处理残基PRO与LEU，如要处理，请仔细查看结果。
可在‘文件输入与参数控制区’调整输入文件名以及需要旋转的二面角以及二面角采样间隔，也可通过 python.py -i -ri -ci -c1 -c2 -c3- c4- -c5 ，后一种可覆盖前一种的参数
注意：采样间隔是度（°）而非弧度（rad）,且都是整数。
