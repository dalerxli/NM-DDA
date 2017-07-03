%此脚本用于生成各文件路径

%main dir
dataRootPath = 'D:\DDAdata\' ;

%sub dir
scaDataPath = [dataRootPath , 'scaData\'] ;
outPutMullerPath = [dataRootPath , 'outPutMuller\'] ;

if isdir( dataRootPath ) == 0
    mkdir( dataRootPath )
end

if isdir( scaDataPath ) == 0
    mkdir( scaDataPath )
end

if isdir( outPutMullerPath ) == 0
    mkdir( outPutMullerPath )
end