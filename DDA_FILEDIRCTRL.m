%�˽ű��������ɸ��ļ�·��

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