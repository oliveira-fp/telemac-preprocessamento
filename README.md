# Pré-processamento Operacional para TELEMAC

![Hydrodynamic Modeling](https://img.shields.io/badge/Hydrodynamic_Modeling-TELEMAC-0055AA?style=for-the-badge)
*Sistema de pré-processamento de dados para modelos hidrodinâmicos TELEMAC 3D*

## 📌 Visão Geral

Este projeto consiste em um conjunto de scripts Fortran para pré-processamento de dados meteorológicos e oceanográficos para alimentação do modelo hidrodinâmico TELEMAC 3D. O sistema realiza:

- 📊 Processamento de dados atmosféricos (ventos do modelo ETA)
- 🌊 Interpolação de dados oceanográficos MOM (SSH, Temperatura, Salinidade, Correntes)
- 🎯 Geração de condições iniciais e de contorno para TELEMAC

## 🛠 Stack Tecnológica

![Fortran](https://img.shields.io/badge/Fortran-90/95-734F96?style=for-the-badge&logo=fortran&logoColor=white)
![Bash](https://img.shields.io/badge/Bash-4EAA25?style=for-the-badge&logo=gnu-bash&logoColor=white)
![NetCDF](https://img.shields.io/badge/NetCDF-3498DB?style=for-the-badge&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4AkEEjEXsR7Z9gAAAB1pVFh0Q29tbWVudAAAAAAAQ3JlYXRlZCB3aXRoIEdJTVBkLmUHAAAAVUlEQVQ4y2NgGAWjFDCqAUa1gQkQfwDi/4zQzJgYqgH+B2IuIP4PxP+Q5KHVQdX8B2IuIP4PxP+Q5KHVQdX8B2IuIP4PxP+Q5KHVQdUMAwB8Fh9ZQm3vQAAAAABJRU5ErkJggg==)
![TELEMAC](https://img.shields.io/badge/TELEMAC-3D-00599C?style=for-the-badge)

**Principais dependências:**
- `gfortran` (compilador Fortran)
- `netcdf-fortran` (biblioteca para manipulação de arquivos NetCDF)
- `cdo` (Climate Data Operators - opcional)
- `Bash` (para automação)

## 📂 Estrutura do Projeto

```bash
scripts/
├── run_prep_grade_4.bash     # Script principal de execução
├── namelist_pre.nml          # Arquivo de configuração
│
├── 10_trata_vento_eta.F90    # Processamento de dados de vento ETA
├── 11_remove_nivel_mom.F90   # Processamento de nível do mar (SSH)
├── 12_cc_sal_temp_mom.F90    # Processamento de salinidade/temperatura
├── 13_cc_velocidade_mom.F90  # Processamento de velocidades
│
└── resultados/               # Diretório de saídas
```

## 🚀 Fluxo de Processamento

### Compilação e Execução
```bash
./n_prep_grade_4.bash  # Compila e executa todos os módulos sequencialmente
```

### Saídas Geradas

**Dados Atmosféricos:**
- `Malha_Vento_ETA03.dat` - Ventos interpolados para a malha

**Condições Iniciais:**
- `Sal_ini.dat`, `Tem_ini.dat` - Salinidade e temperatura iniciais
- `Vel_u_ini.dat`, `Vel_v_ini.dat` - Velocidades iniciais

**Condições de Contorno:**
- `nivel_cont.dat` - Elevação do nível do mar
- `cont_sal_temp.dat` - Salinidade e temperatura
- `cont_velo.dat` - Velocidades

## 🔧 Módulos Principais

### `10_trata_vento_eta.F90` - Processamento de Vento
![Wind Processing](https://img.shields.io/badge/Module-Wind_Processing-00AAFF?style=flat)

- 📥 Leitura de dados do modelo ETA (NetCDF)
- 📍 Interpolação bilinear para pontos da malha TELEMAC
- 💨 Geração de séries temporais de componentes u/v

### `11_remove_nivel_mom.F90` - Processamento de SSH
![SSH Processing](https://img.shields.io/badge/Module-SSH_Processing-0099FF?style=flat)

- 🌊 Extração de SSH (Sea Surface Height) de arquivos MOM NetCDF
- 🎯 Interpolação para pontos de contorno
- ⚡ Tratamento de condições de fronteira aberta

### `12_cc_sal_temp_mom.F90` - Processamento T/S
![TS Processing](https://img.shields.io/badge/Module-T_S_Processing-FF6600?style=flat)

- 🌡️ Interpolação vertical de perfis de temperatura e salinidade
- 🚫 Tratamento de valores NaN/missing (substituição por 1000 → 0)
- 📊 Geração de condições iniciais e de contorno 3D

### `13_cc_velocidade_mom.F90` - Processamento de Correntes
![Currents Processing](https://img.shields.io/badge/Module-Currents_Processing-0066FF?style=flat)

- 🔄 Processamento de componentes u/v da corrente oceânica
- 📏 Interpolação vertical para diferentes níveis sigma
- 🌊 Geração de condições de contorno 3D

## ⚙️ Arquivo de Configuração (`namelist_pre.nml`)

O arquivo de configuração permite personalizar:

```fortran
&atmosferico
  ncfile_eta = "caminho/para/dados/eta.nc"  ! Dados ETA
  xlon = "longitude_variable"               ! Nome da variável de longitude
  ylat = "latitude_variable"                ! Nome da variável de latitude
  uwind = "u10"                             ! Variável de vento U
  vwind = "v10"                             ! Variável de vento V
/

&ssh_config
  ncfile_mom = "caminho/para/dados/mom.nc"  ! Dados MOM
  hsea = "ssh"                              ! Variável de SSH
/

&temp_salt
  ttemp = "temp"                            ! Variável de temperatura
  salt = "salt"                             ! Variável de salinidade
/

&corrente
  uwind = "uo"                              ! Variável de velocidade U
  vwind = "vo"                              ! Variável de velocidade V
/
```

## 💻 Como Executar

### Pré-requisitos
```bash
# Instalar dependências no Ubuntu/Debian
sudo apt-get install gfortran libnetcdff-dev netcdf-bin
```

### Execução Completa
```bash
# 1. Clone o repositório
git clone [URL_DO_REPOSITORIO]
cd preprocessamento

# 2. Configure os caminhos (editar cria_namelist.bash se necessário)
./cria_namelist.bash

# 3. Execute o processamento completo
chmod +x n_prep_grade_4.bash
./n_prep_grade_4.bash
```

### Execução Individual
```bash
# Compilar e executar módulo específico
gfortran 10_trata_vento_eta.F90 -o vento.o -lnetcdff
./vento.o
```

## 📊 Estrutura de Dados

### Arquivos de Entrada Necessários
```
../dados/
├── Eta03_BESM_*.nc          # Dados atmosféricos ETA
├── *ocean_telemac*.nc       # Dados oceanográficos MOM
└── c_contorno_BC.cli        # Definição de contornos

../resultados/
├── malha_2d_deg.dat         # Malha 2D (lat/lon)
├── zz_coluna.dat           # Malha 3D (UTM)
└── zz_coluna_deg.dat       # Malha 3D (lat/lon)
```

### Formatos de Saída
- **Condições Iniciais**: Formato colunar com [ID, X, Y, Valor]
- **Condições de Contorno**: Formato temporal com cabeçalho de nós
- **Precisão**: Float64 para coordenadas, Float32 para dados

## 🎯 Características Técnicas

### Interpolação
- **Horizontal**: Método do vizinho mais próximo
- **Vertical**: Interpolação por profundidade
- **Temporal**: Conservação de séries temporais

### Tratamento de Dados
- ✅ Substituição de valores NaN por 1000 → 0
- ✅ Validação de ranges (UTM, coordenadas)
- ✅ Compressão NetCDF para economizar espaço

### Performance
- ⚡ Processamento eficiente com alocação dinâmica
- 💾 Liberação explícita de memória após uso
- 📦 Processamento por variável para evitar overflow

## 🐛 Solução de Problemas

### Erros Comuns
1. **"Morto" (Out of Memory)**
   ```bash
   # Reduzir uso de memória
   ulimit -s unlimited
   export OMP_NUM_THREADS=1
   ```

2. **Erro NetCDF**
   - Verificar caminhos dos arquivos
   - Confirmar nomes das variáveis no namelist

3. **Problemas de Compilação**
   ```bash
   # Verificar instalação do NetCDF
   nc-config --all
   ```

---

**Desenvolvido para processamento de modelo hidrodinâmico TELEMAC 3D**  

**Autor:** Fernando Pereira de Oliveira  

[![LinkedIn](https://img.shields.io/badge/LinkedIn-0077B5?style=for-the-badge&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/fernando-oliveira-612963245/)
[![GitHub](https://img.shields.io/badge/GitHub-100000?style=for-the-badge&logo=github&logoColor=white)](https://github.com/oliveira-fp)

![Hydrodynamic Modeling](https://img.shields.io/badge/Hydrodynamic_Modeling-TELEMAC_3D-0055AA?style=for-the-badge)



