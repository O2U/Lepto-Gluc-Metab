# Host-like Temperatures Unlock Glucose Metabolism in Pathogenic Leptospira

[![License: MIT](https://img.shields.io/badge/Code%20License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![License: CC BY 4.0](https://img.shields.io/badge/Data%20License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

---

### 概要 (Abstract)

本研究では、病原性細菌レプトスピラが宿主の体温(37°C)をシグナルとして、これまで見過ごされてきたグルコース代謝系を活性化させ、持続感染に利用していることを明らかにしました。ゲノムスケール代謝モデル、¹³C標識プロテオミクス、およびマウス感染モデルを組み合わせることで、この温度依存的な代謝スイッチが本菌の病原性における重要な因子であることを突き止めました。

---

### リポジトリの構成 (Repository Structure)

このリポジトリは、論文の解析を再現するためのコード、データ、結果を格納しています。

- `/code`: 解析に用いた全てのPythonおよびRスクリプト
- `/data`: 生データ、処理済みデータ、および外部データへのリンク
- `/results`: スクリプト実行によって生成される図、表、代謝モデル
- `/manuscript`: 論文のプレプリント

---

### 環境構築 (Setup)

本解析は、Condaを用いた環境管理を推奨します。以下の手順で、必要なソフトウェアとライブラリをインストールできます。

1.  **Condaをインストールします。** ([Miniconda](https://docs.conda.io/en/latest/miniconda.html)を推奨)
2.  **このリポジトリをクローンします。**
    ```bash
    git clone [https://github.com/your-username/Host-Leptospira-Metabolism.git](https://github.com/your-username/Host-Leptospira-Metabolism.git)
    cd Host-Leptospira-Metabolism
    ```
3.  **Conda環境を構築・有効化します。**
    ```bash
    conda env create -f environment.yml
    conda activate leptospira-env
    ```

---

### 解析の再現手順 (Workflow)

以下の順序で`code/`ディレクトリ内のスクリプトを実行することで、論文の主要な結果を再現できます。

1.  **GENRE構築とRIPTiDe解析:**
    ```bash
    cd code/01_genre_construction/
    python run_genre_construction.py
    cd ../02_riptide_analysis/
    python run_riptide.py
    ```
    * **生成物:** `results/models/leptospira_genre.xml`

2.  **In Vitro実験データの解析と作図:**
    ```R
    # RStudioなどで以下のスクリプトを実行
    source("code/03_in_vitro_analysis/analyze_growth_and_qpcr.R")
    ```
    * **生成物:** `results/figures/figure_2.png`

3.  **(以降、プロテオミクス、In Vivo解析、最終的な作図スクリプトの実行方法を同様に記載)**

---

### 大規模データ (Large Data Files)

プロテオミクスのRAWファイルなど、GB単位の大規模データは本リポジトリには含まれていません。全ての生データは、以下のZenodoリポジトリからダウンロード可能です。

- **Zenodo Archive:** [https://doi.org/10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX)

---

### 引用 (Citation)

この研究、コード、またはデータがあなたの研究の役に立った場合は、以下の論文を引用してください。

> (ここに論文の引用情報を記載)

データセットの引用はこちら:

> (ここにZenodoの引用情報を記載)

---

### ライセンス (License)

- **コード:** 本リポジトリ内の全コードは[MITライセンス](LICENSE)の下で公開されています。
- **データと図:** データおよび図は[Creative Commons Attribution 4.0 International (CC-BY 4.0)](https://creativecommons.org/licenses/by/4.0/)ライセンスの下で公開されています。