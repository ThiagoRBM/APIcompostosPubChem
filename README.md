# APIcompostosPubChem
APIcompostosPubChem


Script com função para buscar informações sobre compostos químicos no site PubChem (https://pubchem.ncbi.nlm.nih.gov/) através da própria API deles.
Até agora, são retornadas as seguintes informações: índice de Kovats, número CID, Massa Molecular, Fórmula Molecular e SMILES. É possível adicionar mais coisas.

A função recebe uma coluna com nomes de compostos (em inglês), busca esses nomes no PubChem e pega as informações mencionadas. O reultado é salvo como um arquivo
.txt no diretório indicado, com o nome indicado no início do script. Na tabela salva, tem os compostos encontrados primeiro, com as respectivas informações e os
compostos não encontrados, abaixo.
Caso compostos comuns não sejam encontrados, é provável que eles não estejam escritos de forma correta, ou estejam escritos de alguma forma que o PubChem não reconhece.
