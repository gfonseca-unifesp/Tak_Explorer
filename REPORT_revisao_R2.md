# Relatório — Revisão R2 (ECOLIND-59694R1)

Branch: `revision-r2`. Nenhum push/PR foi feito. Nenhum dado bruto foi alterado
(`ISA_DeepData_2026.csv` e `example_data.csv` continuam byte a byte como estavam).

---

## (i) Tabela 2 verificada — o que pude confirmar e o que não pude

Script de verificação: `scratchpad/analyze_table2*.py` (não versionado, gitignored).
Resultado célula a célula: `tak_table2_verificada.csv` (raiz do repo).

### Achado central

`ISA_DeepData_2026.csv` **não é** uma tabela em nível de registro de ocorrência.
É uma tabela já agregada por linhagem taxonômica: 2.816 linhas, cada uma um
caminho taxonômico único (`Dataset` + Filo..Espécie), com uma única coluna `n`.//
Confirmei que `n` é uma medida de **abundância/densidade** (individuos, ou uma
estimativa de densidade de levantamento fotográfico/de vídeo — não é inteiro em
210 das 2.155 linhas do CCZ e em 58 das 293 da Indian Ocean; é sempre inteiro em
Mid Atlantic e NW Pacific), **não** uma contagem de registros de ocorrência.
Não há nenhuma outra coluna no arquivo que codifique "quantos registros brutos
foram agregados nesta linha".

Consequência prática:

- **A metade "táxons únicos" da Tabela 2 é reconstruível a partir deste arquivo**
  e bate quase exatamente com o publicado: **20 das 24 células batem
  exatamente**; as outras 4 diferem por 1–2 táxons (ver tabela abaixo). É uma
  discrepância pequena, plausivelmente uma sincronização levemente diferente do
  OBIS/DeepData entre a extração usada no artigo e este arquivo — não achei
  nenhuma causa espúria (não é duplicata por maiúscula/espaço: testei
  explicitamente e não há nenhuma).
- **A metade "registros terminais" (os valores entre parênteses) e os totais de
  linha (166.871 / 2.771 / 2.099 / 740) NÃO são reconstruíveis a partir deste
  arquivo.** Testei todas as hipóteses razoáveis — soma de `n` por rank terminal,
  soma de `n` em colunas não-NA, contagem de linhas por rank terminal, razões
  constantes por região — nenhuma reproduz as magnitudes publicadas (as
  contagens de linha ficam ~10×–360× menores que os valores publicados, sem
  fator de escala constante entre ranks). **860 espécies em 22.052 registros**
  (CCZ), por exemplo: os 860 batem (são as 860 linhas terminadas em Espécie),
  mas não existe combinação de colunas neste arquivo que produza 22.052.

  **Conclusão: os "registros terminais" da Tabela 2 do artigo foram calculados
  a partir do pull bruto de ocorrências do OBIS/DeepData (nível de registro,
  não de linhagem agregada) — esse arquivo não está neste repositório.** Se
  quiser que eu valide essa metade da tabela também, preciso do arquivo bruto
  (um registro de ocorrência por linha).

### Reconciliação do total do CCZ (166.871 vs. soma dos parênteses 166.813)

A diferença de 58 **não** pode ser explicada por "registros do CCZ sem
identificação em nenhum rank" *neste* arquivo: há exatamente **1** linha no CCZ
com todos os 6 ranks vazios (n=18, não 58). Como o arquivo não é o pull bruto
de ocorrências (ver acima), esse gap de 58 registros deve estar no arquivo bruto
que gerou a Tabela 2 original, não neste `ISA_DeepData_2026.csv`. Os totais das
4 regiões batem exatamente com o grand total do artigo (166.871+2.771+2.099+740
= 172.481) — isso é consistência interna do próprio artigo, não algo que eu
verifiquei de forma independente a partir deste arquivo.

### Outras afirmações do texto

| Afirmação | Status | Nota |
|---|---|---|
| 36 filos no conjunto | ✅ Confere | 36 valores distintos de `Phylum` no arquivo inteiro |
| CCZ: 860 espécies | ✅ Confere | |
| Indian Ocean: 7 espécies | ✅ Confere | |
| NW Pacific: 3 espécies | ✅ Confere | |
| CCZ: 860 espécies em 22.052 registros | ⚠️ Parcial | 860 confere; 22.052 não é reconstruível (ver acima) |
| 172.481 registros de ocorrência no total | ❌ Não verificável aqui | Soma interna do artigo bate; arquivo não permite conferir de forma independente |
| 2.876 registros taxonômicos únicos | ❌ Não verificável aqui | Testei várias definições candidatas (contagem de linhas=2.816, tuplas de linhagem únicas cross-região=2.476, soma de únicos por rank=3.086) — nenhuma bate com 2.876 |
| East Pacific→NW Pacific, East Indian→Indian Ocean | ✅ Já renomeado | O arquivo já usa os nomes novos (`NW Pacific`, `Indian Ocean`) — nenhum resquício do nome antigo |

### Células "táxons únicos" com pequena divergência (não explicadas por duplicata textual)

| Região | Rank | Publicado | Recalculado |
|---|---|---|---|
| CCZ | Order | 184 | 185 |
| CCZ | Family | 561 | 563 |
| CCZ | Genus | 1094 | 1095 |
| Indian Ocean | Order | 69 | 70 |

Todas as outras 20 células (incluindo todas as de NW Pacific e Mid Atlantic)
batem exatamente.

**Recomendação**: para o manuscrito, cite os números de "táxons únicos" como
confirmados (com a pequena ressalva acima, provavelmente arredondamento de uma
sincronização de dados ligeiramente diferente); para a metade "registros" e os
totais, ou (a) anexe/me envie o pull bruto de ocorrências para eu validar
também, ou (b) mantenha os números como estão (vieram de um cálculo anterior
não reproduzido aqui) e apenas documente no método que vieram do pull bruto,
não do `ISA_DeepData_2026.csv` agregado.

---

## (ii) Explicação do total do CCZ

Ver acima — a diferença de 58 registros não é reconstruível a partir do
`ISA_DeepData_2026.csv` fornecido (esse arquivo só tem 1 linha CCZ sem nenhuma
identificação, não 58). A causa deve estar no pull bruto de ocorrências
original, fora deste repositório.

---

## (iii) Modos de peso: texto do manuscrito × código

O manuscrito diz que o usuário pode escolher pesos "uniform, linear, or
exponential". Isso **não bate integralmente** com nenhuma das duas
implementações:

- **`app.R` (R/Shiny, a ferramenta interativa)**: não existe nenhum seletor de
  modo. Há só um campo de texto livre "Weight Vector" (`textInput`, default
  `"1, 2, 3, 4, 5, 6"`) — o usuário digita os pesos manualmente, número por
  número. "Uniform" (`1,1,1,1,1,1`) e "linear" (`1,2,3,4,5,6`) são alcançáveis
  digitando à mão, mas não existem como opção nomeada/botão; "exponencial" não
  tem fórmula nenhuma implementada — o usuário teria que calcular os valores
  fora do app e colar.
- **`TaK_fun_2.py` (Python, a função standalone)**: tem um parâmetro
  `simple_weight` (booleano) — `True` dá peso uniforme (todos 1), `False` sem
  `custom_weights` dá peso linear (`np.arange(1, n+1)`), e `custom_weights` aceita
  qualquer vetor arbitrário. **"Exponential" não existe em lugar nenhum do
  código** — nem fórmula, nem parâmetro, nem menção em comentário.

**Resumo**: uniform e linear existem como modos nomeados só no lado Python;
exponential não existe em nenhum dos dois. Não implementei uma fórmula
exponencial "chutada" (base 2? `e`? potência do rank?) porque isso é uma
decisão de modelagem dos autores, não uma correção de bug — prefiro que vocês
decidam a fórmula antes de eu codificar. Duas saídas possíveis:

1. **Alinhar o manuscrito ao código**: trocar "uniform, linear, or exponential"
   por algo como "user-defined weights, with uniform and linear as built-in
   presets in the Python reference implementation" (e no R, mencionar que os
   pesos são inteiramente livres via o campo de texto).
2. **Alinhar o código ao manuscrito**: adicionar um modo exponential explícito
   (fórmula e base a definir por vocês) em `TaK_fun_1.R`, `TaK_fun_2.py` e como
   opção no `app.R`. Não fiz essa mudança nesta rodada — avise a fórmula
   desejada e eu implemento.

---

## (iv) Paridade R × Python

**Resultado: paridade confirmada (TR, TC e a tabela-resumo), depois de corrigir
uma divergência real de comportamento.**

Rodei as duas implementações sobre `example_data.csv` e `ISA_DeepData_2026.csv`,
mesmo vetor de pesos `[1,2,3,4,5,6]` nas duas (scripts em
`scratchpad/run_R_parity.R` / `scratchpad/run_py_parity.py`, não versionados).

**Divergência real encontrada e corrigida**: em `ISA_DeepData_2026.csv` há 1
linha do CCZ sem identificação em nenhum rank (Phylum..Species todos vazios).
`dplyr::group_by()` (R) mantém um grupo com `NA` como grupo próprio por padrão
— então o resumo em R tinha 81 linhas. `pandas.groupby()` (Python) **descarta
silenciosamente** linhas com valor `NaN` na coluna de agrupamento por padrão
(`dropna=True`) — o resumo em Python tinha só 80 linhas, perdendo essa linha
inteira sem aviso nenhum. Corrigido em `TaK_fun_2.py` passando
`dropna=False` no `groupby()`, para casar com o comportamento do R.

Depois da correção:

| Dataset | Linhas R | Linhas Python | Máx. diferença TC | Máx. diferença TR | Máx. diferença unique_taxa |
|---|---|---|---|---|---|
| `example_data.csv` | 10 | 10 | 3.3e-16 | 4.4e-16 | 0 |
| `ISA_DeepData_2026.csv` | 81 | 81 | 5.0e-16 | 5.0e-16 | 0 |

As diferenças remanescentes são ruído de ponto flutuante (ordem de 1e-16),
essencialmente idênticas. As fórmulas de TC/TR nas duas linguagens são
matematicamente equivalentes (conferi álgebra: `TC_R = TC_Python` e
`TR_R = TR_Python` termo a termo antes mesmo de rodar) — o único bug real era o
`dropna` do pandas, não a matemática.

---

## (v) O que mudou no app e como testar

### Renomeação de colunas (Tarefa 2.1)

Coluna `n_Lineages` (R) / `n_lines` (Python) → **`unique_taxa`** em
`TaK_fun_1.R`, `app.R` (a função está duplicada inline lá, não faz `source()`
de `TaK_fun_1.R` — atualizei as duas cópias), `Vignette_app.Rmd` e
`TaK_fun_2.py`. O rollup por Dataset `Total_Lineages` → `Total_unique_taxa`
(aba Summary do app e `descriptive_stats` do Python).

**Decisão deliberada — não criei uma coluna `terminal_records`.** A Tarefa 1
mostrou que `Total_N` / `n_individuals` (soma de `Abundance`/`count_col`) **não
é** uma contagem de registros de ocorrência — é abundância/densidade, que pode
inclusive não ser inteira. Chamar essa coluna de "terminal_records" teria
introduzido uma afirmação falsa sobre o que ela mede. Além disso, a
tabela-resumo do app é por `Dataset × primeiro rank taxonômico` (normalmente
Filo) — ela nunca calculou uma quebra por rank ao estilo da Tabela 2
("registros cuja identificação terminal foi em X"); esse é um cálculo
estruturalmente diferente (o que fiz à parte, como script, para a Tarefa 1).
Documentei isso em comentário no topo de `calculate_TaK_shiny()`
(R, nos dois arquivos) e da função `TaK()` (Python).

### Baixa confiança para amostras pequenas (Tarefa 2.4, opcional)

Novo campo na barra lateral do `app.R`: **"Low-confidence threshold (min.
records per group)"** (default 10). Cada grupo `Dataset × rank` com `Total_N`
abaixo do limiar ganha `Confidence = "Low confidence"` (senão `"OK"`) —
refletido no Biplot como círculo vazado (baixa confiança) vs. preenchido (OK).
**O cálculo de TR/TC não muda** — é só sinalização visual, como pedido.

### Bugs reais corrigidos ao longo do caminho

- `Vignette_app.Rmd` fazia `source("TaK_funshiny.R")` — esse arquivo **não
  existe** no repositório (só `TaK_fun_1.R` e `TaK_fun_2.py`). Corrigido para
  `source("TaK_fun_1.R")`. Sem essa correção, renderizar a vignette falha
  imediatamente no primeiro chunk.
- `TaK_fun_2.py`: `groupby()` descartava silenciosamente linhas com
  `group_col` vazio (ver item iv acima) — corrigido com `dropna=False`.

### Como testar

```r
# R: checagem de sintaxe rápida
Rscript -e "invisible(lapply(c('TaK_fun_1.R','app.R'), parse))"

# Rodar o app
shiny::runApp("app.R")
```

1. Aba **Data Editor**: os dados de exemplo embutidos carregam sem erro; subir
   `ISA_DeepData_2026.csv` via "Upload CSV" também carrega sem erro (2.816
   linhas, 4 datasets).
2. Aba **Visualization**: Biplot e "Frequency per Quadrant" renderizam para os
   dois conjuntos de dados; no Biplot, pontos com poucos registros aparecem
   como círculo vazado quando abaixo do limiar configurado na barra lateral.
3. Aba **Summary**: a coluna `Total_unique_taxa` aparece no lugar de
   `Total_Lineages`; para `ISA_DeepData_2026.csv`, `CCZ` mostra
   `Total_unique_taxa = 2155` (bate com o total de linhas do CCZ no arquivo).
4. Download CSV/PNG em cada aba funciona sem erro.

Testei tudo isso de fato rodando o app (R 4.3.2, todas as dependências já
instaladas), incluindo upload real do `ISA_DeepData_2026.csv` via injeção de
arquivo no `<input type=file>` do navegador — sem erros no console do
navegador nem do servidor R em nenhum passo.

### Figuras regeneradas (Tarefa 2.5)

Salvei em `figures/` (novo diretório, versionado):

- `fig1_biplot_example_data.png` — Biplot a partir dos dados de exemplo
  embutidos no app (equivalente à Fig. 1 do artigo).
- `fig2_biplot_ISA_DeepData.png` — Biplot a partir de `ISA_DeepData_2026.csv`
  (equivalente à Fig. 2).
- `fig3_quadrant_ISA_DeepData.png` — Gráfico de frequência por quadrante a
  partir de `ISA_DeepData_2026.csv` (equivalente à Fig. 3).

Geradas pelo próprio app rodando (`downloadHandler` real via `ggsave()`), não
por um script à parte — então refletem exatamente o que o app produz hoje,
incluindo o novo mapeamento de forma por `Confidence`.

### Não fiz nesta rodada

- Não implementei um modo de peso "exponential" (ver item iii — decisão dos
  autores primeiro).
- Não refatorei a duplicação da função `calculate_TaK_shiny()` entre
  `TaK_fun_1.R` e `app.R` (o app não faz `source()` do arquivo, tem uma cópia
  inline) — mantive as duas cópias em sincronia manualmente para esta tarefa,
  mas isso é um risco de manutenção a médio prazo: uma mudança futura em um
  lugar pode não ser replicada no outro. Sugiro considerar `app.R` fazer
  `source("TaK_fun_1.R")` em vez de duplicar a função, como um fast-follow.
