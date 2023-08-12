function Untitled = importfile(filename, dataLines)
%IMPORTFILE1 从文本文件中导入数据
%  UNTITLED = IMPORTFILE1(FILENAME)读取文本文件 FILENAME 中默认选定范围的数据。  以表形式返回数据。
%
%  UNTITLED = IMPORTFILE1(FILE, DATALINES)按指定行间隔读取文本文件 FILENAME
%  中的数据。对于不连续的行间隔，请将 DATALINES 指定为正整数标量或 N×2 正整数标量数组。
%
%  示例:
%  Untitled = importfile1("E:\MyProject\matlab\MATLAB Drive\DtN FEM\Results\case 1&2(FEM and BEM_various materials)\code\Untitled.dat", [1, Inf]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2023-04-26 13:05:59 自动生成

%% 输入处理

% 如果不指定 dataLines，请定义默认范围
if nargin < 2
	dataLines = [1, Inf];
end

%% 设置导入选项并导入数据
opts = delimitedTextImportOptions("NumVariables", 11);

% 指定范围和分隔符
opts.DataLines = dataLines;
opts.Delimiter = ",";

% 指定列名称和类型
opts.VariableNames = ["Type", "Index", "space", "Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8"];
opts.VariableTypes = ["categorical", "double", "string", "double", "double", "double", "double", "double", "double", "double", "double"];

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 指定变量属性
opts = setvaropts(opts, "space", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Type", "space"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["Var4", "Var5", "Var6", "Var7", "Var8"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["Var4", "Var5", "Var6", "Var7", "Var8"], "ThousandsSeparator", ",");

% 导入数据
Untitled = readtable(filename, opts);

end