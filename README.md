# Host-like Temperatures Unlock Glucose Metabolism in Pathogenic Leptospira

[![License: MIT](https://img.shields.io/badge/Code%20License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![License: CC BY 4.0](https://img.shields.io/badge/Data%20License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

---

### 概要 (Abstract)

本研究では、ゲノムスケール代謝モデルを用いた予測から、病原性細菌レプトスピラが宿主の体温(37°C)をシグナルとして、これまで見過ごされてきたグルコース代謝系を活性化させることを明らかにしました。

---

### リポジトリの構成 (Repository Structure)

このリポジトリは、論文の解析を再現するためのコード、データ、結果を格納しています。

- `/code`: 解析に用いた全てのPythonおよびRスクリプト
- `/data`: 生データ、処理済みデータ、および外部データへのリンク
- `/results`: スクリプト実行によって生成される図、表、代謝モデル
- `/manuscript`: 論文のプレプリント

---

### 環境構築 (Setup)

本解析は、Condaを用いた環境管理を推奨します。以下のコマンドで、解析に必要な全てのソフトウェアとライブラリをインストールできます。

~~~bash
conda env create -f environment.yml
~~~

Rの実行環境（バージョン情報など）は`R_session_info.txt`を参照してください。

---

### 引用 (Citation)

この研究、コード、またはデータがあなたの研究の役に立った場合は、以下の論文を引用してください。

> (ここに論文の引用情報を記載)

---

### ライセンス (License)

- **コード:** 本リポジトリ内の全コードは[MITライセンス](LICENSE)の下で公開されています。
- **データと図:** データおよび図は[Creative Commons Attribution 4.0 International (CC-BY 4.0)](https://creativecommons.org/licenses/by/4.0/)ライセンスの下で公開されています。