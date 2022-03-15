### script para buscar informacoes de compostos no PubChem
### 
### 
### 
### 
### 



library(readxl)
library(httr)
library(dplyr)
library(stringr)
library(rlist)

nomesCompostos <- read_excel("gitAPIpubChem/nomesCompostos.xlsx") ## abrir a tabela do exemplo (ou uma outra)

names(nomesCompostos) ### verificar qual a coluna que tem os nomes dos compostos (no caso dessa
## do exemplo que tem só uma coluna, é fácil, mas podem ser tabelas maiores) para os quais se deseja
## informacoes
## 

diretorio= "C:/Users/HP/Desktop/" ## diretorio onde quer que a tabela com as informacoes dos compostos seja salva
## substituir aqui, mantendo as aspas

nomeTabela= "tabelaCompostos2" ## nome que quer que a tabela tenha, no diretorio em que sera criada

###### NAO MEXER AQUI EMBAIXO (INDICO QUANDO FOR PARA MEXER DE NOVO)
###### 
###### 
###### 
###### 
###### 
###### 
###### 
###### 

InfoPubChem= function(tabela, colunaNomesCompostos){ ## roda daqui para baixo sem alterar nada, 
  ## até ser indicado que é para rodar!

print("Iniciando a busca no PubChem")
  
compostos= unique(tabela[[colunaNomesCompostos]]) %>% 
  gsub(" ", "%20", .)
  
url= paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", compostos, "/JSON")

CID= lapply(url, function(x) httr::GET(x, inherits=FALSE)) ## procurando os compostos no PubChem
names(CID)= compostos ## os valores de CID serão usados na busca e não os nomes propriamente ditos
## por causa da maneira como a API foi construida (isso em 03/2022)

print("Buscando as informacoes dos compostos")

NaoEncontradosNome= names(list.filter(CID, status_code != 200))

pagina = lapply(CID, function(x) httr::content(x, encoding = "UTF-8",
                                               type = 'text/plain'))
names(pagina)= compostos

CIDok= unique(sapply(pagina, 
                     function(x) stringr::str_extract(gsub("\n", "", x), '(?<=\"cid\":)\\s*(.*?)\\s*(?=\"atoms\")') %>% 
                       gsub("[^[:digit:].]", "", .))) 
## extraindo o CID do que foi encontrado e pegando só os valores únicos


urlCID= paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/", CIDok, "/JSON")
## montando as URLs para fazer a busca, baseada nos numeros CID

tabCID= lapply(urlCID, function(x) GET(x, inherits=FALSE)) ## acessando o pubChem

tabCID200= list.filter(tabCID, status_code == 200) ## pegando apenas o que tiver dado status 200
## nesses casos, o composto foi encontrado

NaoEncontradosCID= names(list.filter(tabCID, status_code != 200))

tabPlain = lapply(tabCID200, function(x) httr::content(x, encoding = "UTF-8",
                                                          type = 'text/plain'))

print("Montando a tabela com as informacoes")

Lista= vector("list", length= length(unique(tabPlain)))

for (i in 1:length(tabPlain)) {
  
  if (is.null(grepl("Standard non-polar", tabPlain[[i]])) == FALSE){ ## vendo se existe informação
    ## de índice de Kovats com coluna apolar. Se sim, pegar informações
    ## nome do composto, "record numero" (CID, talvez)
    
    
    composto= str_extract(gsub("\n", "",tabPlain[[i]]), "(?<=RecordTitle)\\s*(.*?)\\s*(?=Section)") %>% 
      str_extract(., '(?<=\": \")\\s*(.*?)\\s*(?=\",)')
    
    numero= str_extract(gsub("\n", "",tabPlain[[i]]), "(?<=RecordNumber)\\s*(.*?)\\s*(?=RecordTitle)") %>% 
      gsub("[^[:digit:].]", "", .)
    
    coluna= ifelse(is.na(str_extract(gsub("\n", "",tabPlain[[i]]), "Standard non-polar")),
                   str_extract(gsub("\n", "",tabPlain[[i]]), "Semi-standard non-polar"),
                   str_extract(gsub("\n", "",tabPlain[[i]]), "Standard non-polar"))
    ## se não tiver o rótulo "standard non polar", pegar a coluna "semi standard non polar"
    ## e a partir do nome da coluna, pegar os respectivos kovats (logo abaixo):
    
    ## se o rotulo coluna (gerado acima) nao tiver a palavra "semi", pegar os valores de indice de Kovats
    ## da palavra "standard non polar", em caso de existir a palavra semi, pegar os numeros de kovats da
    ## coluna semi standard non polar
    IK= ifelse(grepl("Semi", coluna, ignore.case=TRUE) == FALSE,
               str_extract(gsub("\n", "", tabPlain[[i]]), "Standard non-polar\\s*(.*?)\\s*(?=ReferenceNumber)") %>%
                 gsub("\\s+", " ", .) %>% 
                 gsub("[^[:digit:].,]", "", .) %>% 
                 gsub(",", " ", .) %>% 
                 str_trim(., side= c("both")) %>% 
                 strsplit(., " ") %>% 
                 as.data.frame(),
               
               str_extract(gsub("\n", "", tabPlain[[i]]), "Semi-standard non-polar\\s*(.*?)\\s*(?=TOCHeading)") %>%
                 gsub("\\s+", " ", .) %>% 
                 gsub("[^[:digit:].,]", "", .) %>% 
                 gsub(",", " ", .) %>% 
                 str_trim(., side= c("both")) %>% 
                 strsplit(., " ") %>% 
                 as.data.frame()
    )
    
    MassaMolecular= str_extract(gsub("\n", "",tabPlain[[i]]), "(?<=Molecular Weight)\\s*(.*?)\\s*(?=g/mol)") %>% 
      str_extract(., "(?<=String)\\s*(.*?)\\s*(?=Unit)") %>% 
      str_extract(., "[0-9.]+") ## pegando a Massa Molecular
    
    Smiles= str_extract(gsub("\n", "",tabPlain[[i]]), "(?<=Canonical SMILES)\\s*(.*?)\\s*(?=\"TOCHeading)") %>% 
      str_extract(., '(?<=String)\\s*(.*?)\\s*(?=\\" +\\})') %>% 
      str_extract(., '(?<=String)\\s*(.*?)\\s*(?=$)') %>% 
      str_extract(., "[CNOH].*") ## pegando os SMILES
    
    FormulaMolecular= str_extract(gsub("\n", "",tabPlain[[i]]), '(Molecular Formula).+?(?=/)') %>% 
      str_extract(., '([C])(\\d+).+?(?=\\")') ## pegando a formula molecular
    
    names(Lista)[i] = composto
    Lista[[i]][["nome"]] = ifelse(is.null(coluna) == TRUE, "NA", composto)
    Lista[[i]][["RecordNumber"]] = ifelse(is.null(coluna) == TRUE, "NA", numero)
    Lista[[i]][["coluna"]] = ifelse(is.null(coluna) == TRUE, "NA", coluna)
    Lista[[i]][["Kovats"]] = ifelse(is.null(IK) == TRUE, "NA", IK)
    Lista[[i]][["MassaMolecular"]] = ifelse(is.null(MassaMolecular) == TRUE, "NA", MassaMolecular)
    Lista[[i]][["SMILES"]] = ifelse(is.null(Smiles) == TRUE, "NA", Smiles)
    Lista[[i]][["FormulaMoecular"]] = ifelse(is.null(FormulaMolecular) == TRUE, "NA", FormulaMolecular)
    
    
  } else { ## se nao tiver coluna non polar, nao pegar o numero de kovats
    
    Lista[[i]] = "Sem kovats"
    
  }
  
  if (i == length(tabPlain)){
    TabTotal= lapply(Lista, function(x) data.frame(name= unlist(x[1]), 
                                                      cid= unlist(x[2]), 
                                                      coluna= unlist(x[3]),
                                                      kovats= unlist(x[4]),
                                                      MassaMolecular=unlist(x[5]),
                                                      FormulaMolecular=unlist(x[7]),
                                                      SMILES=unlist(x[6])
    ))
    
    naoEncontrados= data.frame(name= c(NaoEncontradosCID, NaoEncontradosNome), 
                               cid= "nao encontrado", 
                               coluna= NA,
                               kovats= NA,
                               MassaMolecular= NA,
                               FormulaMolecular= NA,
                               SMILES= NA)
    
    Encontrados= do.call( "rbind", TabTotal) %>% 
      distinct(name, coluna, kovats, .keep_all = TRUE)
    row.names(TabKovats)= NULL
    
    Compostos= rbind(Encontrados, naoEncontrados)
    
  }
  
}

print("Tabela montada")

write.table(Compostos, paste0(diretorio, nomeTabela, ".txt"), sep= "\t",
            fileEncoding = "UTF-16LE", row.names = FALSE, col.names = TRUE)

print(paste0("Tabela salva em: ", diretorio, " como ", nomeTabela, ".txt"))
cat("\n")
print(paste0(diretorio, nomeTabela, ".txt"))

print("Fim")
return(TabKovats)


}


##### MEXER ABAIXO:

infoCompostos= InfoPubChem(nomesCompostos, "Composto") ## substituir "nomesCompostos" pelo nome que dado à tabela carregada
## no R (SEM AS ASPAS) e substituir "Composto" pelo nome da coluna em que os nomes dos compostos para os quais
## informações devem ser buscadas (COM AS ASPAS)
