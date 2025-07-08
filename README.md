# Pré-processamento Operacional para TELEMAC

![Hydrodynamic Modeling](https://img.shields.io/badge/Hydrodynamic_Modeling-TELEMAC-0055AA?style=for-the-badge)
*Sistema de pré-processamento de dados para modelos hidrodinâmicos*

## 📌 Visão Geral

Este projeto consiste em um conjunto de scripts para pré-processamento de dados meteorológicos e oceanográficos para alimentação do modelo hidrodinâmico TELEMAC. O sistema realiza:

- Processamento de dados atmosféricos (ventos)
- Interpolação de dados oceanográficos (SSH, Temperatura, Salinidade, Correntes)
- Geração de condições iniciais e de contorno
- Conversão entre sistemas de coordenadas

## 🛠 Stack Tecnológica

![Fortran](https://img.shields.io/badge/Fortran-%23734F96.svg?style=for-the-badge&logo=fortran&logoColor=white)
![Bash](https://img.shields.io/badge/Bash-4EAA25?style=for-the-badge&logo=gnu-bash&logoColor=white)
![NetCDF](https://img.shields.io/badge/NetCDF-3498DB?style=for-the-badge&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4AkEEjEXsR7Z9gAAAB1pVFh0Q29tbWVudAAAAAAAQ3JlYXRlZCB3aXRoIEdJTVBkLmUHAAAAVUlEQVQ4y2NgGAWjFDCqAUa1gQkQfwDi/4zQzJgYqgH+B2IuIP4PxP+Q5KHVQdX8B2IuIP4PxP+Q5KHVQdX8B2IuIP4PxP+Q5KHVQdX8B2IuIP4PxP+Q5KHVQdUMAwB8Fh9ZQm3vQAAAAABJRU5ErkJggg==)
![TELEMAC](https://img.shields.io/badge/TELEMAC-00599C?style=for-the-badge&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4AkEEjEXsR7Z9gAAAB1pVFh0Q29tbWVudAAAAAAAQ3JlYXRlZCB3aXRoIEdJTVBkLmUHAAABJElEQVQ4y2NgGAWjFDCqAUa1gQkQfwDi/4zQzJgYqgH+B2IuIP4PxP+Q5KHVQdX8B2IuIP4PxP+Q5KHVQdX8B2IuIP4PxP+Q5KHVQdX8B2IuIP4PxP+Q5KHVQdUMAwB8Fh9ZQm3vQAAAAABJRU5ErkJggg==)

**Principais dependências:**
- `gfortran` (compilador Fortran)
- `netcdf-fortran` (biblioteca para manipulação de arquivos NetCDF)
- Bash (para automação)

## 📂 Estrutura do Projeto

```bash
preprocessamento/
├── run_pre.bash            # Script principal de execução
├── cria_namelist.bash      # Gerador de arquivo de configuração
├── 1_trata_vento_era.F90   # Processamento de dados de vento
├── 2_remove_nivel_mom.F90  # Processamento de nível do mar
├── 3_cc_sal_temp_mom.F90   # Processamento de salinidade/temperatura
├── 4_cc_velocidade_mom.F90 # Processamento de velocidades
└── namelist_pre.nml        # Arquivo de configuração gerado
```

## 🚀 Fluxo de Processamento

1. **Configuração Inicial**
   ```bash
   ./cria_namelist.bash  # Gera arquivo de configuração
   ```

2. **Compilação e Execução**
   ```bash
   ./run_pre.bash  # Compila e executa todos os módulos
   ```

3. **Saídas Geradas**
   - `Malha_Vento_ETA03.dat` (dados de vento interpolados)
   - `nivel_cont.dat` (condições de contorno para SSH)
   - `Sal_ini.dat`, `Tem_ini.dat` (condições iniciais T/S)
   - `Vel_u_ini.dat`, `Vel_v_ini.dat` (condições iniciais de corrente)
   - `cont_sal_temp.dat`, `cont_velo.dat` (condições de contorno)

## 🔧 Módulos Principais

### 1. Processamento de Vento (`1_trata_vento_era.F90`)
![Wind Processing](https://img.shields.io/badge/Module-Wind_Processing-00AAFF?style=flat)

- Interpolação bilinear de campos de vento
- Geração de séries temporais para pontos da malha

### 2. Processamento de Nível do Mar (`2_remove_nivel_mom.F90`)
![SSH Processing](https://img.shields.io/badge/Module-SSH_Processing-0099FF?style=flat)

- Extração de SSH (Sea Surface Height) de arquivos NetCDF
- Interpolação para pontos de contorno
- Tratamento de condições de fronteira aberta

### 3. Processamento T/S (`3_cc_sal_temp_mom.F90`)
![TS Processing](https://img.shields.io/badge/Module-T_S_Processing-FF6600?style=flat)

- Interpolação vertical de perfis de temperatura e salinidade
- Tratamento de valores NaN/missing
- Geração de condições iniciais e de contorno

### 4. Processamento de Correntes (`4_cc_velocidade_mom.F90`)
![Currents Processing](https://img.shields.io/badge/Module-Currents_Processing-0066FF?style=flat)

- Processamento de componentes u/v da corrente
- Interpolação para diferentes níveis de profundidade
- Geração de condições de contorno para fronteiras abertas

## 💻 Como Executar

1. Clone o repositório:
```bash
git clone [URL_DO_REPOSITORIO]
cd preprocessamento
```

2. Configure os caminhos dos arquivos de entrada no `cria_namelist.bash`

3. Execute o processamento completo:
```bash
chmod +x run_pre.bash
./run_pre.bash
```

## 📊 Dependências de Dados

**Arquivos de entrada necessários:**
- Dados atmosféricos (NetCDF): `../dados/Eta03_BESM_*.nc`
- Dados oceanográficos (NetCDF): `../dados/*ocean_telemac*.nc`
- Malha computacional: `../resultados/malha_deg.dat`
- Definição de contornos: `../dados/c_contorno_BC.cli`

## 🛠️ Personalização

Edite o arquivo `namelist_pre.nml` para ajustar:
- Nomes de variáveis nos arquivos NetCDF
- Parâmetros temporais
- Configurações específicas do domínio

## 📝 Licença

Este projeto está licenciado sob a Licença MIT - veja o arquivo [LICENSE](LICENSE) para detalhes.

---
Desenvolvido como parte do processamento para modelos hidrodinâmicos: TELEMAC.  

**Autor:** Fernando Pereira de Oliveira  

[![LinkedIn](https://img.shields.io/badge/LinkedIn-0077B5?style=for-the-badge&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/fernando-oliveira-612963245/)
[![GitHub](https://img.shields.io/badge/GitHub-100000?style=for-the-badge&logo=github&logoColor=white)](https://github.com/oliveira-fp)

![Hydrodynamic Modeling](https://img.shields.io/badge/Hydrodynamic_Modeling-TELEMAC-0055AA?style=for-the-badge)
